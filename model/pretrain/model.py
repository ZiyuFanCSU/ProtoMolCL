import numpy as np
from rdkit import DataStructs
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch_geometric.utils import to_dense_batch
from torch.nn.utils.rnn import pad_sequence
from sklearn.metrics.pairwise import cosine_similarity
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem
from model.models.aggr import PromptAggr

class LayerSelector(nn.Module):
    def __init__(self, channels, num_layers):
        super().__init__()
        self.gate = nn.Linear(channels, 1)
        self.num_layers = num_layers

    def forward(self, outs):
        h = torch.stack(outs, dim=0)
        score = self.gate(h)
        score = score.squeeze(-1).transpose(0, 1)
        weights = torch.softmax(score, dim=-1)
        selected = torch.sum(h.transpose(0, 1) * weights.unsqueeze(-1), dim=1)
        return selected, weights

class GNNWrapper(nn.Module):
    def __init__(self, gnn, emb_dim, config, atom_context_size, aggr_head=4, layer_norm=True, layer_norm_out=True):
        super(GNNWrapper, self).__init__()
        self.device = config['device']
        self.margin = config['optim']['margin']
        self.dist_metric = config['optim']['distance_metric']
        self.num_candidates = config['optim']['num_candidates']
        self.dropout_ratio = config['model']['dropout_ratio']
        self.temperature = 0.1
        self.INF = 1e8
        self.ssl_batch = config['batch_size_ssl']
        self.prompt_token = ['<low>', '<mid>', '<high>']
        self.prompt_inds = torch.LongTensor(
            [0] + [1] * 5  + [2] * (self.num_candidates - 1) + [3])
        self.emb_dim = emb_dim
        self.gnn = gnn
        self.aggrs = nn.ModuleList([PromptAggr(emb_dim=emb_dim,
                                               num_heads=aggr_head,
                                               dropout=self.dropout_ratio,
                                               local_loss=i==(len(self.prompt_token)-1),
                                               layer_norm_out=layer_norm_out)
                                    for i in range(len(self.prompt_token))])
        self.context_func = nn.Sequential(
            nn.Linear(emb_dim, emb_dim), nn.ReLU(),
            nn.Linear(emb_dim, atom_context_size)
        )
        self.motif_func = nn.Sequential(
            nn.Linear(emb_dim, emb_dim), nn.ReLU(),
            nn.Linear(emb_dim, emb_dim), nn.ReLU(),
            nn.Linear(emb_dim, 86)
        )
        self.proj_low_inter = nn.Sequential(
            nn.Linear(emb_dim, emb_dim),
            nn.ReLU(),
            nn.Linear(emb_dim, emb_dim)
        )
        self.proj_mid_inter = nn.Sequential(
            nn.Linear(emb_dim, emb_dim),
            nn.ReLU(),
            nn.Linear(emb_dim, emb_dim)
        )
        self.proj_high_inter = nn.Sequential(
            nn.Linear(emb_dim, emb_dim),
            nn.ReLU(),
            nn.Linear(emb_dim, emb_dim)
        )

        self.proj_low_intra = nn.Sequential(
            nn.Linear(emb_dim, emb_dim),
            nn.ReLU(),
            nn.Linear(emb_dim, emb_dim)
        )

        self.proj_mid_intra = nn.Sequential(
            nn.Linear(emb_dim, emb_dim),
            nn.ReLU(),
            nn.Linear(emb_dim, emb_dim)
        )

        self.proj_high_intra = nn.Sequential(
            nn.Linear(emb_dim, emb_dim),
            nn.ReLU(),
            nn.Linear(emb_dim, emb_dim)
        )

        self.selector_low = LayerSelector(emb_dim, 2)
        self.selector_mid = LayerSelector(emb_dim, 3)
        self.selector_high = LayerSelector(emb_dim, 3)

        self.proto_low = nn.Sequential(
            nn.Linear(emb_dim, emb_dim),
            nn.ReLU(),
            nn.Linear(emb_dim, emb_dim)
        )

        self.proto_mid = nn.Sequential(
            nn.Linear(emb_dim, emb_dim),
            nn.ReLU(),
            nn.Linear(emb_dim, emb_dim)
        )

        self.proto_high = nn.Sequential(
            nn.Linear(emb_dim, emb_dim),
            nn.ReLU(),
            nn.Linear(emb_dim, emb_dim)
        )

        self.prototypes_low = nn.Parameter(torch.randn(34, emb_dim))
        self.prototypes_mid = nn.Parameter(torch.randn(108, emb_dim))
        self.prototypes_high = nn.Parameter(torch.randn(8, emb_dim))

        nn.init.xavier_uniform_(self.prototypes_low)
        nn.init.xavier_uniform_(self.prototypes_mid)
        nn.init.xavier_uniform_(self.prototypes_high)

    def prototype_similarity(self, h, prototypes, tau=0.2):
        h = F.normalize(h, dim=-1)
        prototypes = F.normalize(prototypes, dim=-1)
        logits = h @ prototypes.t()
        probs = torch.softmax(logits / tau, dim=-1)
        return logits, probs
            
    def distance(self, tensor_a, tensor_b, metric=None):
        if (metric is None and self.dist_metric == 'cossim') or metric == 'cossim':
            return 1 - F.cosine_similarity(tensor_a, tensor_b, dim=-1)
        elif (metric is None and self.dist_metric == 'l2norm') or metric == 'l2norm':
            return (tensor_a - tensor_b).norm(dim=-1)
        else:
            raise Exception
        
    def get_graphs_low(self,graph_reps,batch):
        batch_x_0, batch_mask_0 = to_dense_batch(graph_reps["re_low"][0], batch.batch)
        batch_x_1, batch_mask_1 = to_dense_batch(graph_reps["re_low"][1], batch.batch)
        batch_size = batch_x_0.size(0) // len(self.prompt_inds)
        prompt_inds = torch.concat([self.prompt_inds] * batch_size).to(batch_x_0.device)
        graph_lows_0=self.aggrs[0](batch_x_0[prompt_inds == 0], batch_mask_0[prompt_inds == 0])
        graph_lows_1 = self.aggrs[0](batch_x_1[prompt_inds == 0], batch_mask_1[prompt_inds == 0])
        out_low = [graph_lows_0[0],graph_lows_1[0]]
        g_low, w_low = self.selector_low(out_low)
        return g_low, w_low, graph_lows_0, graph_lows_1
    
    def get_graphs_mid(self,graph_reps,batch):
        batch_x_0, batch_mask_0 = to_dense_batch(graph_reps["re_mid"][0], batch.batch)
        batch_x_1, batch_mask_1 = to_dense_batch(graph_reps["re_mid"][1], batch.batch)
        batch_x_2, batch_mask_2 = to_dense_batch(graph_reps["re_mid"][2], batch.batch)
        batch_size = batch_x_0.size(0) // len(self.prompt_inds)
        prompt_inds = torch.concat([self.prompt_inds] * batch_size).to(batch_x_0.device)
        graph_mids_0=self.aggrs[1](batch_x_0[prompt_inds == 0], batch_mask_0[prompt_inds == 0])
        graph_mids_1 = self.aggrs[1](batch_x_1[prompt_inds == 0], batch_mask_1[prompt_inds == 0])
        graph_mids_2 = self.aggrs[1](batch_x_2[prompt_inds == 0], batch_mask_2[prompt_inds == 0])
        out_mid = [graph_mids_0[0],graph_mids_1[0],graph_mids_2[0]]
        g_mid, w_mid = self.selector_mid(out_mid)
        return g_mid, w_mid, graph_mids_0, graph_mids_1, graph_mids_2
    
    def get_graphs_high(self,graph_reps,batch):
        batch_x_0, batch_mask_0 = to_dense_batch(graph_reps["re_high"][0], batch.batch)
        batch_x_1, batch_mask_1 = to_dense_batch(graph_reps["re_high"][1], batch.batch)
        batch_x_2, batch_mask_2 = to_dense_batch(graph_reps["re_high"][2], batch.batch)
        batch_size = batch_x_0.size(0) // len(self.prompt_inds)
        prompt_inds = torch.concat([self.prompt_inds] * batch_size).to(batch_x_0.device)
        graph_highs_0=self.aggrs[2](batch_x_0[prompt_inds == 0], batch_mask_0[prompt_inds == 0])
        graph_highs_1 = self.aggrs[2](batch_x_1[prompt_inds == 0], batch_mask_1[prompt_inds == 0])
        graph_highs_2 = self.aggrs[2](batch_x_2[prompt_inds == 0], batch_mask_2[prompt_inds == 0])

        out_high = [graph_highs_0[0],graph_highs_1[0],graph_highs_2[0]]
        g_high, w_high = self.selector_high(out_high)
        return g_high, w_high, graph_highs_0, graph_highs_1, graph_highs_2
    
    def get_graphs_high_mask(self,graph_reps,batch):
        batch_x_0, batch_mask_0 = to_dense_batch(graph_reps["re_high"][0], batch.batch)
        batch_x_1, batch_mask_1 = to_dense_batch(graph_reps["re_high"][1], batch.batch)
        batch_x_2, batch_mask_2 = to_dense_batch(graph_reps["re_high"][2], batch.batch)
        batch_x_0=batch_x_0[:self.ssl_batch*len(self.prompt_inds)]
        batch_mask_0=batch_mask_0[:self.ssl_batch*len(self.prompt_inds)]
        batch_x_1=batch_x_1[:self.ssl_batch*len(self.prompt_inds)]
        batch_mask_1=batch_mask_1[:self.ssl_batch*len(self.prompt_inds)]
        batch_x_2=batch_x_2[:self.ssl_batch*len(self.prompt_inds)]
        batch_mask_2=batch_mask_2[:self.ssl_batch*len(self.prompt_inds)]
        batch_size = self.ssl_batch
        prompt_inds = torch.concat([self.prompt_inds] * batch_size).to(batch_x_0.device)
        graph_highs_0=self.aggrs[2](batch_x_0[prompt_inds == 3], batch_mask_0[prompt_inds == 3])
        graph_highs_1 = self.aggrs[2](batch_x_1[prompt_inds == 3], batch_mask_1[prompt_inds == 3])
        graph_highs_2 = self.aggrs[2](batch_x_2[prompt_inds == 3], batch_mask_2[prompt_inds == 3])
        out_high = [graph_highs_0[0],graph_highs_1[0],graph_highs_2[0]]
        g_high, w_high = self.selector_high(out_high)
        return g_high, w_high, graph_highs_0, graph_highs_1, graph_highs_2
    
    def get_graphs_low_inter(self,graph_reps,batch):
        batch_x_0, batch_mask_0 = to_dense_batch(graph_reps["re_low"][0], batch.batch)
        batch_x_1, batch_mask_1 = to_dense_batch(graph_reps["re_low"][1], batch.batch)
        batch_x_0=batch_x_0[:self.ssl_batch*len(self.prompt_inds)]
        batch_mask_0=batch_mask_0[:self.ssl_batch*len(self.prompt_inds)]
        batch_x_1=batch_x_1[:self.ssl_batch*len(self.prompt_inds)]
        batch_mask_1=batch_mask_1[:self.ssl_batch*len(self.prompt_inds)]
        batch_size = self.ssl_batch
        prompt_inds = torch.concat([self.prompt_inds] * batch_size).to(batch_x_0.device)
        graph_lows_0=self.aggrs[0](batch_x_0[(prompt_inds == 0) | (prompt_inds == 2)], batch_mask_0[(prompt_inds == 0) | (prompt_inds == 2)])
        graph_lows_1 = self.aggrs[0](batch_x_1[(prompt_inds == 0) | (prompt_inds == 2)], batch_mask_1[(prompt_inds == 0) | (prompt_inds == 2)])
        out_low = [graph_lows_0[0],graph_lows_1[0]]
        g_low, w_low = self.selector_low(out_low)
        return g_low, w_low, graph_lows_0, graph_lows_1
    
    def get_graphs_mid_inter(self,graph_reps,batch):
        batch_x_0, batch_mask_0 = to_dense_batch(graph_reps["re_mid"][0], batch.batch)
        batch_x_1, batch_mask_1 = to_dense_batch(graph_reps["re_mid"][1], batch.batch)
        batch_x_2, batch_mask_2 = to_dense_batch(graph_reps["re_mid"][2], batch.batch)
        batch_x_0=batch_x_0[:self.ssl_batch*len(self.prompt_inds)]
        batch_mask_0=batch_mask_0[:self.ssl_batch*len(self.prompt_inds)]
        batch_x_1=batch_x_1[:self.ssl_batch*len(self.prompt_inds)]
        batch_mask_1=batch_mask_1[:self.ssl_batch*len(self.prompt_inds)]
        batch_x_2=batch_x_2[:self.ssl_batch*len(self.prompt_inds)]
        batch_mask_2=batch_mask_2[:self.ssl_batch*len(self.prompt_inds)]
        batch_size = self.ssl_batch
        prompt_inds = torch.concat([self.prompt_inds] * batch_size).to(batch_x_0.device)
        graph_mids_0=self.aggrs[1](batch_x_0[(prompt_inds == 0) | (prompt_inds == 2)], batch_mask_0[(prompt_inds == 0) | (prompt_inds == 2)])
        graph_mids_1 = self.aggrs[1](batch_x_1[(prompt_inds == 0) | (prompt_inds == 2)], batch_mask_1[(prompt_inds == 0) | (prompt_inds == 2)])
        graph_mids_2 = self.aggrs[1](batch_x_2[(prompt_inds == 0) | (prompt_inds == 2)], batch_mask_2[(prompt_inds == 0) | (prompt_inds == 2)])
        out_mid = [graph_mids_0[0],graph_mids_1[0],graph_mids_2[0]]
        g_mid, w_mid = self.selector_mid(out_mid)
        return g_mid, w_mid, graph_mids_0, graph_mids_1, graph_mids_2
    
    def get_graphs_high_inter(self,graph_reps,batch):
        batch_x_0, batch_mask_0 = to_dense_batch(graph_reps["re_high"][0], batch.batch)
        batch_x_1, batch_mask_1 = to_dense_batch(graph_reps["re_high"][1], batch.batch)
        batch_x_2, batch_mask_2 = to_dense_batch(graph_reps["re_high"][2], batch.batch)
        batch_x_0=batch_x_0[:self.ssl_batch*len(self.prompt_inds)]
        batch_mask_0=batch_mask_0[:self.ssl_batch*len(self.prompt_inds)]
        batch_x_1=batch_x_1[:self.ssl_batch*len(self.prompt_inds)]
        batch_mask_1=batch_mask_1[:self.ssl_batch*len(self.prompt_inds)]
        batch_x_2=batch_x_2[:self.ssl_batch*len(self.prompt_inds)]
        batch_mask_2=batch_mask_2[:self.ssl_batch*len(self.prompt_inds)]
        batch_size = self.ssl_batch
        prompt_inds = torch.concat([self.prompt_inds] * batch_size).to(batch_x_0.device)
        graph_highs_0=self.aggrs[2](batch_x_0[(prompt_inds == 0) | (prompt_inds == 2)], batch_mask_0[(prompt_inds == 0) | (prompt_inds == 2)])
        graph_highs_1 = self.aggrs[2](batch_x_1[(prompt_inds == 0) | (prompt_inds == 2)], batch_mask_1[(prompt_inds == 0) | (prompt_inds == 2)])
        graph_highs_2 = self.aggrs[2](batch_x_2[(prompt_inds == 0) | (prompt_inds == 2)], batch_mask_2[(prompt_inds == 0) | (prompt_inds == 2)])
        out_high = [graph_highs_0[0],graph_highs_1[0],graph_highs_2[0]]
        g_high, w_high = self.selector_high(out_high)
        return g_high, w_high, graph_highs_0, graph_highs_1, graph_highs_2

    def get_representation(self, batch):
        graph_reps = self.gnn(batch.x, batch.pe, batch.edge_index, batch.edge_attr, batch.batch)
        g_low, w_low, graph_lows_0,graph_lows_1 = self.get_graphs_low(graph_reps,batch)
        g_mid, w_mid, graph_mids_0, graph_mids_1,graph_mids_2= self.get_graphs_mid(graph_reps,batch)
        g_high, w_high, graph_highs_0, graph_highs_1, graph_highs_2 = self.get_graphs_high(graph_reps,batch)

        return g_low, w_low, graph_lows_0,graph_lows_1, g_mid, w_mid, graph_mids_0, graph_mids_1,graph_mids_2,g_high, w_high, graph_highs_0, graph_highs_1, graph_highs_2, graph_reps

    def forward(self, batch):
        g_low, w_low, graph_lows_0,graph_lows_1, g_mid, w_mid, graph_mids_0, graph_mids_1,graph_mids_2,g_high, w_high, graph_highs_0, graph_highs_1, graph_highs_2, graph_reps = self.get_representation(batch)
        proto_loss, proto_loss_dict = self.compute_loss_proto(g_low, w_low, g_mid, w_mid, g_high, w_high, batch)
        numhighlayer = round(proto_loss_dict["depth_high"]-3)
        g_high_mask, w_high_mask, graph_highs_0_mask, graph_highs_1_mask, graph_highs_2_mask = self.get_graphs_high_mask(graph_reps,batch)
        if numhighlayer == 0:
            loss_highmask = self.compute_loss_highmask(batch, graph_highs_0_mask)
        elif numhighlayer == 1: 
            loss_highmask = self.compute_loss_highmask(batch, graph_highs_1_mask)
        elif numhighlayer == 2: 
            loss_highmask = self.compute_loss_highmask(batch, graph_highs_2_mask)
        loss_lowmotif = self.compute_loss_lowmotif(batch,g_low)
        loss_2, loss_1_reg_mid_0, loss_1_reg_mid_1, loss_1_reg_mid_2= self.compute_loss_intralayer_mid(graph_reps, batch)
        loss_1_reg_mid = loss_1_reg_mid_0+ loss_1_reg_mid_1+ loss_1_reg_mid_2
        loss_level_align,loss_lm,loss_mh,loss_lh,acc_lm,acc_mh,acc_lh = self.compute_loss_interlayer(graph_reps, batch)
        loss = {"proto_loss":proto_loss,
                "loss_proto_low": proto_loss_dict["loss_proto_low"],
                "loss_proto_mid": proto_loss_dict["loss_proto_mid"],
                "loss_proto_high": proto_loss_dict["loss_proto_high"],
                "depth_low":proto_loss_dict["depth_low"],
                "depth_mid":proto_loss_dict["depth_mid"],
                "depth_high":proto_loss_dict["depth_high"],
                "mask_loss":loss_highmask,
                "motif_loss":loss_lowmotif,
                "margin_loss_mid":loss_2,
                "loss_1_reg_mid":loss_1_reg_mid ,
                "layerinter_loss":loss_level_align,
                "layerinter_loss_lm":loss_lm.detach().item(),
                "layerinter_loss_mh":loss_mh.detach().item(),
                "layerinter_loss_lh":loss_lh.detach().item(),
                "acc_lm":acc_lm.detach().item(),
                "acc_mh":acc_mh.detach().item(),
                "acc_lh":acc_lh.detach().item(),
                }
        return loss

    def compute_loss_interlayer(self, graph_reps, batch):
        g_low, w_low, _,_ = self.get_graphs_low_inter(graph_reps,batch)
        g_mid, w_mid, _,_,_ = self.get_graphs_mid_inter(graph_reps,batch)
        g_high, w_high, _,_,_ = self.get_graphs_high_inter(graph_reps,batch)
        z_low = self.proj_low_inter(g_low)
        z_mid = self.proj_mid_inter(g_mid)
        z_high = self.proj_high_inter(g_high)
        loss_lm, acc_lm = self.infonce_one_side(z_low, z_mid)
        loss_ml, acc_ml = self.infonce_one_side(z_mid, z_low)
        loss_lm = (loss_lm + loss_ml) / 2
        acc_lm = (acc_lm + acc_ml) / 2
        loss_mh, acc_mh = self.infonce_one_side(z_mid, z_high)
        loss_hm, acc_hm = self.infonce_one_side(z_high, z_mid)
        loss_mh = (loss_mh + loss_hm) / 2
        acc_mh = (acc_mh + acc_hm) / 2
        loss_lh, acc_lh = self.infonce_one_side(z_low, z_high)
        loss_hl, acc_hl = self.infonce_one_side(z_high, z_low)
        loss_lh = (loss_lh + loss_hl) / 2
        acc_lh = (acc_lh + acc_hl) / 2
        loss_level_align = loss_lm + loss_mh + loss_lh
        return loss_level_align,loss_lm,loss_mh,loss_lh,acc_lm,acc_mh,acc_lh

    def get_score(self,graph_lows,batch_size):
        aggr_score = graph_lows[2]
        aggr_score = aggr_score.view(batch_size, 6, -1)
        aggr_score_regu = aggr_score.clone().detach()
        aggr_score_regu[:] = (1 / (aggr_score > 0).sum(dim=-1, keepdim=True))
        aggr_score_regu[aggr_score == 0] = 0
        loss_1_reg = F.smooth_l1_loss(
            aggr_score.reshape(-1, aggr_score.size(-1)),
            aggr_score_regu.reshape(-1, aggr_score_regu.size(-1)))
        return loss_1_reg

    def compute_loss_intralayer_mid(self,graph_reps, batch):
        batch_x_0, batch_mask_0 = to_dense_batch(graph_reps["re_mid"][0], batch.batch)
        batch_x_1, batch_mask_1 = to_dense_batch(graph_reps["re_mid"][1], batch.batch)
        batch_x_2, batch_mask_2 = to_dense_batch(graph_reps["re_mid"][2], batch.batch)
        batch_x_0=batch_x_0[:self.ssl_batch*len(self.prompt_inds)]
        batch_mask_0=batch_mask_0[:self.ssl_batch*len(self.prompt_inds)]
        batch_x_1=batch_x_1[:self.ssl_batch*len(self.prompt_inds)]
        batch_mask_1=batch_mask_1[:self.ssl_batch*len(self.prompt_inds)]
        batch_x_2=batch_x_2[:self.ssl_batch*len(self.prompt_inds)]
        batch_mask_2=batch_mask_2[:self.ssl_batch*len(self.prompt_inds)]
        batch_size = self.ssl_batch
        prompt_inds = torch.concat([self.prompt_inds] * batch_size).to(batch_x_0.device)
        graph_mids_0=self.aggrs[1](batch_x_0[(prompt_inds == 0) | (prompt_inds == 1)], batch_mask_0[(prompt_inds == 0) | (prompt_inds == 1)])
        graph_mids_1 = self.aggrs[1](batch_x_1[(prompt_inds == 0) | (prompt_inds == 1)], batch_mask_1[(prompt_inds == 0) | (prompt_inds == 1)])
        graph_mids_2 = self.aggrs[1](batch_x_2[(prompt_inds == 0) | (prompt_inds == 1)], batch_mask_2[(prompt_inds == 0) | (prompt_inds == 1)])
        out_mid = [graph_mids_0[0],graph_mids_1[0],graph_mids_2[0]]
        g_mid, w_mid = self.selector_mid(out_mid)
        mol_fps = [batch.mol_fps[i][0] for i in range(len(batch.mol_fps)) if i % (len(self.prompt_inds)) == 0][:self.ssl_batch]
        mol_dist = 1 - torch.Tensor(
            [DataStructs.BulkTanimotoSimilarity(mol_fps[i], mol_fps) for i in range(len(mol_fps))]).to(batch_x_0.device)
        h_gr = self.proj_mid_intra(g_mid)
        h_gr = h_gr.view(batch_size, 6, self.emb_dim)
        loss_1 = self.adaptive_margin_loss(F.normalize(h_gr, dim=-1), mol_dist)
        loss_1_reg_0 = self.get_score(graph_mids_0,batch_size)
        loss_1_reg_1 = self.get_score(graph_mids_1,batch_size)
        loss_1_reg_2 = self.get_score(graph_mids_2,batch_size)
        return loss_1, loss_1_reg_0, loss_1_reg_1, loss_1_reg_2
    
        
    def compute_loss_highmask(self, batch,graph_rep):
        context_labels = [torch.LongTensor(batch.context_label[i]) for i in range(len(batch.context_label)) if
                          i % (len(self.prompt_inds)) == len(self.prompt_inds) - 1][:self.ssl_batch]
        context_labels = pad_sequence(context_labels, batch_first=True, padding_value=0).to(graph_rep[0].device)
        batch_size = self.ssl_batch
        h_tk = graph_rep[1]
        h_tk = h_tk.view(batch_size, 1, -1, self.emb_dim)
        h_tk = h_tk[:, -1]
        h_tk = h_tk[:, :context_labels.size(1)].contiguous().view(-1, self.emb_dim)
        context_logits = self.context_func(h_tk)
        context_loss = F.cross_entropy(context_logits, context_labels.view(-1), reduction='none')
        context_loss = (context_loss[context_labels.view(-1) != 0]).mean()
        return context_loss

    def compute_loss_lowmotif(self, batch,graph_rep):
        motif_labels = torch.FloatTensor(
            np.array([np.array(batch.mol_mds[i][0][-86:]) for i in range(len(batch.mol_mds)) if
             i % (len(self.prompt_inds)) == 0])).to(
            graph_rep.device)[:self.ssl_batch]
        h_gr = graph_rep[:self.ssl_batch]
        motif_scores = self.motif_func(h_gr)
        motif_loss = F.smooth_l1_loss(motif_scores, motif_labels.float())
        return motif_loss
    
    def infonce(self, z1, z2):
        N = len(z1)
        sim_zii= (z1 @ z1.T) / self.temperature 
        sim_zjj = (z2 @ z2.T) / self.temperature
        sim_zij = (z1 @ z2.T) / self.temperature
        sim_zii = sim_zii - self.INF * torch.eye(N, device=z1.device)
        sim_zjj = sim_zjj - self.INF * torch.eye(N, device=z1.device)
        sim_Z = torch.cat([
            torch.cat([sim_zij, sim_zii], dim=1),
            torch.cat([sim_zjj, sim_zij.T], dim=1)], dim=0)
        log_sim_Z = F.log_softmax(sim_Z, dim=1)
        loss = - torch.diag(log_sim_Z).mean()
        with torch.no_grad():
            pred = torch.argmax(sim_zij, dim=1)
            correct = pred.eq(torch.arange(N, device=z1.device)).sum()
            acc = 100 * correct / N
        return loss, acc
    
    def infonce_one_side(self,z1, z2, temperature=0.2):
        z1 = F.normalize(z1, dim=-1)
        z2 = F.normalize(z2, dim=-1)

        N = z1.size(0)
        logits = z1 @ z2.t() / temperature
        labels = torch.arange(N, device=z1.device)

        loss = F.cross_entropy(logits, labels)

        with torch.no_grad():
            pred = torch.argmax(logits, dim=1)
            acc = 100 * pred.eq(labels).float().mean()

        return loss, acc
    
    def infonce_one_side_weight(self,z1, z2, temperature=0.2):
        z1 = F.normalize(z1, dim=-1)
        z2 = F.normalize(z2, dim=-1)

        N = z1.size(0)
        logits = z1 @ z2.t() / temperature
        labels = torch.arange(N, device=z1.device)

        loss = F.cross_entropy(logits, labels)

        with torch.no_grad():
            pred = torch.argmax(logits, dim=1)
            acc = 100 * pred.eq(labels).float().mean()

        return loss, acc
    
    def soft_label_kl_loss(self, pred_prob, target_prob, eps=1e-8):
        target_prob = target_prob.float().to(pred_prob.device)
        target_prob = target_prob / (target_prob.sum(dim=-1, keepdim=True) + eps)
        log_pred = torch.log(pred_prob + eps)
        loss = F.kl_div(
            log_pred,
            target_prob,
            reduction="batchmean"
        )

        return loss
    def prototype_diversity_loss(self, prototypes):
        p = F.normalize(prototypes, dim=-1)
        sim = p @ p.t()
        K = prototypes.size(0)
        eye = torch.eye(K, device=prototypes.device)
        loss = ((sim - eye) ** 2).mean()
        return loss
    
    def selector_depth_expectation(self, w, layer_ids):
        layer_ids = torch.tensor(layer_ids, device=w.device).float()
        depth = (w * layer_ids).sum(dim=-1)  # [B]
        return depth
    
    def depth_order_loss(self, depth_low, depth_mid, depth_high, margin=0.5):
        loss_lm = F.relu(depth_low + margin - depth_mid).mean()
        loss_mh = F.relu(depth_mid + margin - depth_high).mean()
        loss_lh = F.relu(depth_low + margin - depth_high).mean()
        return loss_lm + loss_mh + loss_lh

    def compute_loss_proto(self, g_low, w_low, g_mid, w_mid, g_high, w_high, batch):
        pg_low = self.proto_low(g_low)
        pg_mid = self.proto_mid(g_mid)
        pg_high = self.proto_high(g_high)
        logits_low, prob_low = self.prototype_similarity(
            pg_low, self.prototypes_low
        )
        logits_mid, prob_mid = self.prototype_similarity(
            pg_mid, self.prototypes_mid
        )
        logits_high, prob_high = self.prototype_similarity(
            pg_high, self.prototypes_high
        )

        target_low = torch.tensor(batch.layer1[::len(self.prompt_inds)]).float().to(g_low.device)
        target_mid = torch.tensor(batch.layer3[::len(self.prompt_inds)]).float().to(g_low.device)
        target_high = torch.tensor(batch.layer5[::len(self.prompt_inds)]).float().to(g_low.device)
        loss_proto_low = self.soft_label_kl_loss(prob_low, target_low)
        loss_proto_mid = self.soft_label_kl_loss(prob_mid, target_mid)
        loss_proto_high = self.soft_label_kl_loss(prob_high, target_high)
        loss_proto_target = loss_proto_low + loss_proto_mid + loss_proto_high
        loss_proto_div = (
            self.prototype_diversity_loss(self.prototypes_low) +
            self.prototype_diversity_loss(self.prototypes_mid) +
            self.prototype_diversity_loss(self.prototypes_high)
        )
        depth_low = self.selector_depth_expectation(w_low, [1, 2])
        depth_mid = self.selector_depth_expectation(w_mid, [2, 3, 4])
        depth_high = self.selector_depth_expectation(w_high, [3, 4, 5])
        loss_selector_order = self.depth_order_loss(depth_low, depth_mid, depth_high)

        total_loss = (
            1.0 * loss_proto_target
            + 1 * loss_proto_div
            + 0.3 * loss_selector_order
        )
        loss_dict = {
            "loss_total": total_loss.detach().item(),
            "loss_proto_target": loss_proto_target.detach().item(),
            "loss_proto_low": loss_proto_low.detach().item(),
            "loss_proto_mid": loss_proto_mid.detach().item(),
            "loss_proto_high": loss_proto_high.detach().item(),
            "loss_proto_div": loss_proto_div.detach().item(),
            "loss_selector_order": loss_selector_order.detach().item(),
            "depth_low": depth_low.detach().mean().item(),
            "depth_mid": depth_mid.detach().mean().item(),
            "depth_high": depth_high.detach().mean().item(),
        }

        return total_loss, loss_dict

    def margin_loss(self, graph_rep, margin_factor=None, anchor_inds=None):
        num_candidates = 6
        batch_size = graph_rep.size(0)
        if anchor_inds is not None:
            sample_mask = torch.zeros((batch_size, num_candidates)).bool()
            sample_mask[(torch.arange(batch_size), anchor_inds)] = True
            anchor = graph_rep[sample_mask].unsqueeze(1)
            pos = graph_rep[~sample_mask].view(batch_size, num_candidates - 1, -1)
            graph_rep = torch.concat([anchor, pos], dim=1)

        anchor = graph_rep[:, :1]
        pos = graph_rep[:, 1:]
        pos_sim = self.distance(anchor, pos)
        pos_sim = pos_sim.repeat_interleave(batch_size - 1, dim=0).view(batch_size, batch_size - 1,
                                                                        -1) 

        ori_rep = graph_rep[:, 0]
        interleave = ori_rep.repeat_interleave(batch_size, dim=0).view(batch_size, batch_size, -1)
        repeat = ori_rep.repeat(batch_size, 1).view(batch_size, batch_size, -1)
        sim_matrix = self.distance(interleave, repeat)
        diag_idx = torch.eye(batch_size, dtype=torch.bool)
        sim_matrix = sim_matrix[~diag_idx].view(batch_size, batch_size - 1) 
        neg_sim = sim_matrix.repeat_interleave(num_candidates - 1, dim=-1).view(batch_size, batch_size - 1,
                                                                                num_candidates - 1)
        if margin_factor is not None:
            diag_idx = torch.eye(batch_size, dtype=torch.bool)
            margin_factor = margin_factor[~diag_idx].view(batch_size, batch_size - 1)
            margin_factor_edge = torch.ones((batch_size, batch_size - 1, num_candidates - 1),
                                            device=margin_factor.device)
            margin_factor_edge[margin_factor == 0] = 0
            margin_factor = margin_factor.repeat_interleave(num_candidates - 1, dim=-1).view(batch_size, batch_size - 1,
                                                                                             num_candidates - 1)
            margin = margin_factor * self.margin
        else:
            margin = self.margin
            margin_factor_edge = torch.ones((batch_size, batch_size - 1, 5),
                                            device=graph_rep.device)

        margin_loss = torch.maximum(torch.tensor(0).to(self.device),
                                    margin_factor_edge * (margin + pos_sim - neg_sim)).mean()

        return margin_loss

    def adaptive_margin_loss(self, graph_rep, scaff_dist, anchor_inds=None):
        batch_size = graph_rep.size(0)
        num_candidates = graph_rep.size(1)
        if anchor_inds is not None:
            sample_mask = torch.zeros((batch_size, num_candidates)).bool()
            sample_mask[(torch.arange(batch_size), anchor_inds)] = True
            anchor = graph_rep[sample_mask].unsqueeze(1)
            pos = graph_rep[~sample_mask].view(batch_size, num_candidates - 1, -1)
            graph_rep = torch.concat([anchor, pos], dim=1)

        margin_loss = self.margin_loss(graph_rep, margin_factor=scaff_dist)
        batch_size = graph_rep.size(0)
        scaff_dist_diff = scaff_dist.unsqueeze(1) - scaff_dist.unsqueeze(2) 

        inds_triplet = torch.argwhere(scaff_dist_diff > 0.3)
        inds_triplet = inds_triplet[(inds_triplet[:, 0] != inds_triplet[:, 1]) &
                                    (inds_triplet[:, 0] != inds_triplet[:, 2])]

        diff_values = scaff_dist_diff[inds_triplet[:, 0], inds_triplet[:, 1], inds_triplet[:, 2]]
        if len(diff_values):
            if len(diff_values) > batch_size ** 2:
                sub_inds = torch.argsort(diff_values, descending=True)[:batch_size ** 2]
                diff_values = diff_values[sub_inds]
                inds_triplet = inds_triplet[sub_inds]

            anchor_inds, pos_inds, neg_inds = inds_triplet[:, 0], inds_triplet[:, 1], inds_triplet[:, 2]
            anchor_rep, pos_rep, neg_rep = graph_rep[anchor_inds, 0], graph_rep[pos_inds, 0], graph_rep[neg_inds, 0]

            pos_sim = self.distance(anchor_rep, pos_rep)
            neg_sim = self.distance(anchor_rep, neg_rep)

            margin_loss_adaptive = torch.maximum(torch.tensor(0).to(graph_rep.device),
                                                 diff_values.detach() * self.margin + pos_sim - neg_sim).mean()

            return margin_loss_adaptive + margin_loss
        else:
            return margin_loss

