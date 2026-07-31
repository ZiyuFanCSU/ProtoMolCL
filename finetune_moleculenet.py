import os
import copy
import random
import numpy as np
import yaml
from tqdm import tqdm
import torch
import torch.nn as nn
from torch_geometric.loader import DataLoader
from torch.utils.data import Subset
from sklearn.metrics import roc_auc_score, r2_score
from torch.optim.lr_scheduler import CosineAnnealingLR
from model.finetune.loader import MoleculeDataset
from model.finetune.model import GNNPredictor
from model.splitters import scaffold_split, moleculeace_split
from model.utils.scheduler import PolynomialDecayLR
from itertools import product

def get_optimizer(model, lr_params):
    assert isinstance(lr_params, dict)
    pretrain_name, finetune_name = [], []
    for name, param in model.named_parameters():
        if 'gnn' in name or 'aggr' in name:
            pretrain_name.append(name)
        elif 'graph_pred_linear' in name:
            finetune_name.append(name)
        else:
            pretrain_name.append(name)

    pretrain_params = list(
        map(lambda x: x[1], list(filter(lambda kv: kv[0] in pretrain_name, model.named_parameters()))))
    finetune_params = list(
        map(lambda x: x[1], list(filter(lambda kv: kv[0] in finetune_name, model.named_parameters()))))

    optimizer = torch.optim.Adam([
        {'params': finetune_params},
        {'params': pretrain_params, 'lr': float(lr_params['pretrain_lr'])}
    ], lr=float(lr_params['finetune_lr']), weight_decay=float(lr_params['decay']))

    return optimizer


def get_dataloader(config, seed=0):
    # Setup dataset
    dataset = MoleculeDataset(config['dataset']['data_dir'],
                              config['dataset']['data_name'],
                              config['dataset']['feat_type'])

    num_task = dataset.num_task
    print('Loading dataset {} of size {} with num_task={}'.format(config['dataset']['data_name'], len(dataset), num_task))

    if 'CHEMBL' in config['dataset']['data_name']:  # MoleculeACE stratified random split
        train_idx, val_idx, test_idx = moleculeace_split(dataset.smiles, dataset.labels, val_size=0.1, test_size=0.1)
    else:  # MoleculeNet scaffold split
        train_idx, val_idx, test_idx = scaffold_split(dataset.smiles, frac_valid=0.1, frac_test=0.1, balanced=False, seed=42)
    
    train_dataset, val_dataset, test_dataset = \
        Subset(dataset, train_idx), Subset(dataset, val_idx), Subset(dataset, test_idx)

    train_loader = DataLoader(train_dataset, batch_size=config['batch_size'], shuffle=True)
    val_loader = DataLoader(val_dataset, batch_size=config['batch_size'], shuffle=False)
    test_loader = DataLoader(test_dataset, batch_size=config['batch_size'], shuffle=False)

    return dataset, train_loader, val_loader, test_loader


def set_seed(seed):
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)


def eval(model, val_loader, config, metric='r2'):
    assert metric in ['rmse', 'r2']
    model.eval()
    y_true, y_scores = [], []
    for step, batch in enumerate(val_loader):
        batch = batch.to(config['device'])
        with torch.no_grad():
            predict = model(batch)['predict']

        y_true.append(batch.label.view(predict.shape))
        y_scores.append(predict)

    y_true = torch.cat(y_true, dim=0).cpu().numpy()
    y_scores = torch.cat(y_scores, dim=0).cpu().numpy()

    if 'CHEMBL' in config['dataset']['data_name']:
        score = r2_score(y_true, y_scores)
    else:
        roc_list = []
        for i in range(y_true.shape[1]):
            # AUC is only defined when there is at least one positive data.
            if np.sum(y_true[:, i] == 1) > 0 and np.sum(y_true[:, i] == -1) > 0:
                is_valid = y_true[:, i] ** 2 > 0
                roc_list.append(roc_auc_score((y_true[is_valid, i] + 1) / 2, y_scores[is_valid, i]))

        score = np.mean(roc_list)

    return score


def train(model, train_loader, criterion, optimizer, scheduler, config, channel_idx=-1):
    model.train()
    loss_history = []
    channel_weight = 0
    for idx, batch in enumerate(train_loader):
        batch.to(config['device'])
        output = model(batch, channel_idx=channel_idx)
        predict = output['predict']
        label = batch.label.view(predict.shape)

        if isinstance(criterion, nn.BCEWithLogitsLoss):
            mask = label == 0  # nan entry
            loss = criterion(predict.double(), (label + 1) / 2) * (~mask)
            loss = loss.sum() / (~mask).sum()
        elif isinstance(criterion, nn.MSELoss):
            loss = criterion(predict, label)
            loss = loss.mean()
        else:
            raise Exception

        optimizer.zero_grad()
        loss.backward()
        if config['optim']['gradient_clip'] > 0:
            nn.utils.clip_grad_norm_(model.parameters(), config['optim']['gradient_clip'])
        optimizer.step()

        if config['optim']['scheduler'] == 'poly_decay':
            scheduler.step()

        loss_history.append(loss.item())

    channel_weight = channel_weight / len(train_loader)

    return np.mean(loss_history), channel_weight




def main(config):
    runseeds = np.random.randint(100, size=config['num_run'])

    # Setup model
    if config['dataset']['feat_type'] == 'basic':
        atom_feat_dim, bond_feat_dim = None, None
    elif config['dataset']['feat_type'] == 'rich':
        atom_feat_dim, bond_feat_dim = 143, 14
    elif config['dataset']['feat_type'] == 'super_rich':
        atom_feat_dim, bond_feat_dim = 170, 14
    else:
        raise NotImplementedError('Unrecognized feature type. Please choose from [basic/rich/super_rich].')


    dataset, train_loader, val_loader, test_loader = get_dataloader(config)
    avg_auc_last, avg_auc_best = [], []

    for i in range(config['num_run']):
        # Setup model
        model = GNNPredictor(num_layer=config['model']['num_layer'],
                             num_layer_low=config['model']['num_layer_low'],
                             num_layer_mid=config['model']['num_layer_mid'],
                             num_layer_high=config['model']['num_layer_high'],
                             emb_dim=config['model']['emb_dim'],
                             num_tasks=dataset.num_task,
                             atom_feat_dim=atom_feat_dim,
                             bond_feat_dim=bond_feat_dim,
                             drop_ratio=config['model']['dropout_ratio'],
                             attn_drop_ratio=config['model']['attn_dropout_ratio'],
                             model_head=config['model']['heads'],
                             layer_norm_out=config['model']['layernorm'])

        if config['model']['checkpoint']:
            checkpoint = torch.load(config['model']['checkpoint'])
            ckpt_state = checkpoint['wrapper']

            model_state = model.state_dict()

            missing_keys, unexpected_keys, shape_mismatch = [], [], []

            for k, v in model_state.items():
                if k in ckpt_state:
                    if ckpt_state[k].shape != v.shape:
                        shape_mismatch.append(k)
                else:
                    missing_keys.append(k)

            for k in ckpt_state.keys():
                if k not in model_state:
                    unexpected_keys.append(k)

            print("====")
            print(config['model']['checkpoint'])
            print('\n'.join(missing_keys))
            print('Loading checkpoint from {}'.format(config['model']['checkpoint']))

            model.load_state_dict(ckpt_state, strict=False)
        model.to(config['device'])


        # Setup optimizer
        optimizer = get_optimizer(model, config['optim'])
        scheduler = None
        if config['optim']['scheduler'] == 'cos_anneal':
            scheduler = CosineAnnealingLR(optimizer, T_max=config['epochs'], eta_min=0.0001)
        elif config['optim']['scheduler'] == 'poly_decay':
            scheduler = PolynomialDecayLR(optimizer, warmup_updates=config['epochs'] * len(train_loader) // 10,
                                          tot_updates=config['epochs'] * len(train_loader),
                                          lr=config['optim']['finetune_lr'], end_lr=1e-9, power=1)

        # Setup loss function
        if config['dataset']['task'] == 'regression':
            criterion = nn.MSELoss(reduction='none')
        elif config['dataset']['task'] == 'classification':
            criterion = nn.BCEWithLogitsLoss(reduction="none")
        else:
            raise NotImplementedError

        best_score, best_checkpoint = -float('inf'), None

        # Setup random seed
        print("Seed:", runseeds[i])
        set_seed(runseeds[i])

        for epoch in tqdm(range(1, config['epochs'] + 1)):
            # train one epoch
            train(model, train_loader, criterion, optimizer, scheduler, config)

            # evaluate validation
            score = eval(model, val_loader, config)
            test_score = eval(model, test_loader, config)

            if config['optim']['scheduler'] == 'cos_anneal':
                scheduler.step()

            cur_finetune_lr = optimizer.param_groups[0]['lr']
            cur_pretrain_lr = optimizer.param_groups[1]['lr']

            tqdm.write(
                f"[ep{epoch}] val={score:>4.4f} test={test_score:>4.4f} "
                f"finetune_lr={cur_finetune_lr:.6g} "
                f"pretrain_lr={cur_pretrain_lr:.6g} "
            )

            if score > best_score:
                best_score = score
                best_checkpoint = copy.deepcopy(model.state_dict())

        score_last_checkpoint = eval(model, test_loader, config)
        avg_auc_last.append(score_last_checkpoint)
        model.load_state_dict(best_checkpoint)
        score_best_checkpoint = eval(model, test_loader, config)
        avg_auc_best.append(score_best_checkpoint)
        
        if 'CHEMBL' in config['dataset']['data_name']:
            print('[Best R2]: {:.4f} {:.4f} {:.4f}'.format(best_score, score_last_checkpoint, score_best_checkpoint))
        else:
            print('[Best AUC]: {:.4f} {:.4f} {:.4f}'.format(best_score, score_last_checkpoint, score_best_checkpoint))

    print(avg_auc_last)
    print('[Last] {} {}'.format(np.mean(avg_auc_last), np.std(avg_auc_last)))
    print(avg_auc_best)
    print('[Best] {} {}'.format(np.mean(avg_auc_best), np.std(avg_auc_best)))

    save_path = config.get('save_path', './result.txt')

    with open(save_path, 'w') as f:
        for i in range(len(avg_auc_best)):
            f.write('Run #{} (seed={}): best={} last={}\n'.format(
                i + 1, runseeds[i], avg_auc_best[i], avg_auc_last[i]
            ))
        f.write('Average last score: {}\n'.format(np.mean(avg_auc_last)))
        f.write('Std last score: {}\n'.format(np.std(avg_auc_last)))
        f.write('Average best score: {}\n'.format(np.mean(avg_auc_best)))
        f.write('Std best score: {}\n'.format(np.std(avg_auc_best)))

    return {
        "last_mean": np.mean(avg_auc_last),
        "last_std": np.std(avg_auc_last),
        "best_mean": np.mean(avg_auc_best),
        "best_std": np.std(avg_auc_best),
    }


if __name__ == '__main__':
    import os
    import copy
    import csv

    tasks=['clintox','bbbp','sider','toxcast','tox21','bace','muv','hiv']
    batch_size_s = [32,64,128]
    epoch_s = [15,20,30,60,100]
    dropout_ratio_s = [0,0.1,0.3]
    attn_dropout_ratio_s = [0,0.1,0.3]
    pretrain_lr_s = [0.0001,0.0005]
    finetune_lr_s = [0.001,0.0005,0.00003]
    batch_size_s = [16]
    result_root = '/result/moleculeace'
    summary_csv = os.path.join(result_root, 'grid_search_summary.csv')

    with open(summary_csv, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow([
            'task',
            'batch_size',
            'epochs',
            'dropout_ratio',
            'attn_dropout_ratio',
            'pretrain_lr',
            'finetune_lr',
            'best_mean',
            'best_std',
            'last_mean',
            'last_std',
            'save_path'
        ])

        for name in tasks:
            yaml_path = './config/moleculenet/' + str(name) + '.yaml'
            with open(yaml_path, 'r') as f:
                base_config = yaml.load(f, Loader=yaml.FullLoader)
            param_grid = product(
                batch_size_s,
                epoch_s,
                dropout_ratio_s,
                attn_dropout_ratio_s,
                pretrain_lr_s,
                finetune_lr_s
            )

            for batch_size, epochs, dropout_ratio, attn_dropout_ratio, pretrain_lr, finetune_lr in param_grid:
                config = copy.deepcopy(base_config)
                config['batch_size'] = batch_size
                config['epochs'] = epochs
                config['dataset']['data_name'] = name
                config['dataset']['feat_type'] = 'super_rich'
                config['model']['checkpoint'] = "./checkpoint/zinc-gps.pt"
                config['model']['dropout_ratio'] = dropout_ratio
                config['model']['attn_dropout_ratio'] = attn_dropout_ratio
                config['optim']['pretrain_lr'] = pretrain_lr
                config['optim']['finetune_lr'] = finetune_lr

                exp_name = (
                    f"{name}"
                    f"_bs{batch_size}"
                    f"_ep{epochs}"
                    f"_drop{dropout_ratio}"
                    f"_attndrop{attn_dropout_ratio}"
                    f"_plr{pretrain_lr}"
                    f"_flr{finetune_lr}"
                )

                task_dir = os.path.join(result_root, "4_" + name)
                os.makedirs(task_dir, exist_ok=True)
                save_path = os.path.join(task_dir, exp_name + '.txt')
                config['save_path'] = save_path

                print("=" * 80)
                print("Current experiment:", exp_name)
                print("=" * 80)

                try:
                    result = main(config)
                    writer.writerow([
                        name,
                        batch_size,
                        epochs,
                        dropout_ratio,
                        attn_dropout_ratio,
                        pretrain_lr,
                        finetune_lr,
                        result['best_mean'],
                        result['best_std'],
                        result['last_mean'],
                        result['last_std'],
                        save_path
                    ])
                    csvfile.flush()

                except Exception as e:
                    print(f"Experiment failed: {exp_name}")
                    print(e)

                    writer.writerow([
                        name,
                        batch_size,
                        epochs,
                        dropout_ratio,
                        attn_dropout_ratio,
                        pretrain_lr,
                        finetune_lr,
                        'ERROR',
                        'ERROR',
                        'ERROR',
                        'ERROR',
                        save_path
                    ])
                    csvfile.flush()

