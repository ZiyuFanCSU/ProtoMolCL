import yaml
import torch
import torch.optim as optim
from torch.utils.data import DataLoader
from model.pretrain.loader import MoleculeDataset
from model.pretrain.model import GNNWrapper
from model.models.gnn import GPS_multilayer
from model.utils.data import list_collate_fn
import warnings
from rdkit import RDLogger
RDLogger.DisableLog("rdApp.warning")
warnings.filterwarnings("ignore", category=DeprecationWarning)

def train_one_epoch(epoch, model, loader, optimizer, config, enable_kg=True):
    model.train()

    accum_iter = config['optim']['accum_iter']
    total_loss_meter = 0.0
    
    pbar = enumerate(loader)

    for batch_idx, batch in pbar:
        batch = batch.to("cuda:0")
        with torch.set_grad_enabled(True):
            loss_dict = model(batch)
            loss = loss_dict["proto_loss"]+loss_dict["mask_loss"]+loss_dict["motif_loss"]+\
                loss_dict["margin_loss_mid"]+config['optim']['reg_coeff'] * loss_dict["loss_1_reg_mid"]+\
                loss_dict["layerinter_loss"]
                     
            loss_reg = loss_dict["loss_1_reg_mid"]
                
            loss = loss / accum_iter
            loss.backward()

            if ((batch_idx + 1) % accum_iter == 0) or (batch_idx + 1 == len(loader)):
                optimizer.step()
                optimizer.zero_grad()
            
            total_loss_meter += loss.item()

            if batch_idx % 20 == 0:
                print(
                    f"[Epoch {epoch}] "
                    f"[Batch {batch_idx}/{len(loader)}] "
                    f"proto=({loss_dict['loss_proto_low']:.2f} | "
                    f"{loss_dict['loss_proto_mid']:.2f} | "
                    f"{loss_dict['loss_proto_high']:.2f}) | "
                    f"margin=({loss_dict['margin_loss_mid'].detach().item():.2f}) | "
                    f"inter=({loss_dict['layerinter_loss_lm']:.2f} | "
                    f"{loss_dict['layerinter_loss_mh']:.2f} | "
                    f"{loss_dict['layerinter_loss_lh']:.2f}) | "
                    f"mask={loss_dict['mask_loss'].detach().item():.3f} | "
                    f"motif={loss_dict['motif_loss'].detach().item():.3f} | "
                    f"loss_reg={loss_reg.detach().item():.2f} | "
                    f"acc=({loss_dict['acc_lm']:.2f} | "
                    f"{loss_dict['acc_mh']:.2f} | "
                    f"{loss_dict['acc_lh']:.2f}) | "
                    f"depth=({loss_dict['depth_low']:.2f}, "
                    f"{loss_dict['depth_mid']:.2f}, "
                    f"{loss_dict['depth_high']:.2f})"
                )

    avg_loss = total_loss_meter / len(loader)
    return avg_loss


def main(config):
    dataset = MoleculeDataset(config['dataset']['data_dir'],
                              num_candidates=config['optim']['num_candidates'],
                              feat_type=config['dataset']['feat_type'])

    loader = DataLoader(dataset, batch_size=config['batch_size'],
                        num_workers=config['dataset']['num_workers'],collate_fn=list_collate_fn,
                        shuffle=False, drop_last=True)
    
    atom_feat_dim, bond_feat_dim = 170, 14

    
    gnn = GPS_multilayer(channels=config['model']['emb_dim'], pe_dim=20,
              node_dim=atom_feat_dim,
              edge_dim=bond_feat_dim,
              num_layers_low=config['model']['num_layer_low'],
              num_layers_mid=config['model']['num_layer_mid'],
              num_layers_high=config['model']['num_layer_high'],
              heads=config['model']['heads'], 
              attn_type='multihead',
              dropout=config['model']['dropout_ratio'],
              attn_dropout=config['model']['dropout_ratio'])

    model = GNNWrapper(gnn, config['model']['emb_dim'], config,
                       atom_context_size=len(dataset.atom_vocab_itos),
                       layer_norm_out=True)

    model = model.to("cuda:0")

    if config['model']['checkpoint']:
        print('Loading wrapper checkpoint from {}'.format(config['model']['checkpoint']))
        loaded_state_dict = torch.load(config['model']['checkpoint'])['wrapper']
        model.load_state_dict(loaded_state_dict, strict=True)

    optimizer = optim.Adam(model.parameters(), lr=config['optim']['lr'], weight_decay=float(config['optim']['decay']))
    
    pbar = range(config['start_epoch'] + 1, config['epochs'] + 1)

    for epoch in pbar:
        loss_total = train_one_epoch(epoch, model, loader, optimizer, config, enable_kg=epoch>20)

        checkpoint_path = '{}/zinc-{}-{}.pt'.format(config['output_dir'],config['model']['backbone'], epoch)
        torch.save({'wrapper': model.state_dict()}, checkpoint_path)
        print('Save to {} at epoch {}'.format(checkpoint_path, epoch))


if __name__ == '__main__':
    with open('./config/pretrain.yaml', 'r') as f:
        config = yaml.load(f, Loader=yaml.FullLoader)
    main(config)