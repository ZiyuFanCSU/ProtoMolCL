import math
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch_geometric.utils import to_dense_batch
from model.models.gnn import GPS_multilayer
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

class GNNPredictor(nn.Module):
    def __init__(self, num_layer,num_layer_low, num_layer_mid,num_layer_high,emb_dim, num_tasks, JK="last", drop_ratio=0, 
                 attn_drop_ratio=0, model_head=4, aggr_head=4,
                 atom_feat_dim=None, bond_feat_dim=None, baseline=None, add_mean_pool=False,
                 layer_norm=True, layer_norm_out=True, act='softmax'):
        super(GNNPredictor, self).__init__()
        self.num_layer = num_layer
        self.num_layer_low = num_layer_low
        self.num_layer_mid = num_layer_mid
        self.num_layer_high = num_layer_high
        self.drop_ratio = drop_ratio
        self.attn_drop_ratio = attn_drop_ratio
        self.JK = JK
        self.emb_dim = emb_dim
        self.num_tasks = num_tasks
        self.baseline = baseline
        self.add_mean_pool = add_mean_pool
        self.layer_norm = layer_norm
        self.layer_norm_out = layer_norm_out
        self.act = act

        if self.num_layer < 2:
            raise ValueError("Number of GNN layers must be greater than 1.")

        
        self.gnn = GPS_multilayer(channels=self.emb_dim, pe_dim=20,
                       node_dim=atom_feat_dim, edge_dim=bond_feat_dim,
                       num_layers_low=self.num_layer_low,
                       num_layers_mid=self.num_layer_mid,
                       num_layers_high=self.num_layer_high, heads=model_head,
                       attn_type='multihead', dropout=self.drop_ratio,
                       attn_dropout=self.attn_drop_ratio)

        self.graph_pred_linear = nn.Sequential(
            nn.Linear(self.emb_dim, self.emb_dim),
            nn.ReLU(),
            nn.Linear(self.emb_dim, self.num_tasks)
        )
        self.attn = nn.Sequential(
            nn.Linear(self.emb_dim, self.emb_dim // 2),
            nn.Tanh(),
            nn.Linear(self.emb_dim // 2, 1)
        )
        self.pred_head = nn.Sequential(
            nn.LayerNorm(self.emb_dim),
            nn.Dropout(self.drop_ratio),
            nn.Linear(self.emb_dim, self.emb_dim),
            nn.ReLU(),
            nn.Dropout(self.drop_ratio),
            nn.Linear(self.emb_dim, num_tasks)
        )

        self.dropout = nn.Dropout(self.drop_ratio)
        self.selector_low = LayerSelector(emb_dim, 2)
        self.selector_mid = LayerSelector(emb_dim, 3)
        self.selector_high = LayerSelector(emb_dim, 3)
        self.dropout_ratio = 0

        self.prompt_token = ['<low>', '<mid>', '<high>']
        self.aggrs = nn.ModuleList([PromptAggr(emb_dim=emb_dim,
                                               num_heads=aggr_head,
                                               dropout=self.dropout_ratio,
                                               local_loss=i==(len(self.prompt_token)-1),
                                               layer_norm_out=layer_norm_out)
                                    for i in range(len(self.prompt_token))])
    
    def freeze_aggr_module(self):
        for param in self.aggrs.parameters():
            param.requires_grad = False

    def get_graphs_low(self,graph_reps,batch):
        batch_x_0, batch_mask_0 = to_dense_batch(graph_reps["re_low"][0], batch.batch)
        batch_x_1, batch_mask_1 = to_dense_batch(graph_reps["re_low"][1], batch.batch)

        graph_lows_0=self.aggrs[0](batch_x_0, batch_mask_0)
        graph_lows_1 = self.aggrs[0](batch_x_1, batch_mask_1)

        out_low = [graph_lows_0[0],graph_lows_1[0]]
        g_low, w_low = self.selector_low(out_low)
        return g_low, w_low, graph_lows_0, graph_lows_1
    
    def get_graphs_mid(self,graph_reps,batch):
        batch_x_0, batch_mask_0 = to_dense_batch(graph_reps["re_mid"][0], batch.batch)
        batch_x_1, batch_mask_1 = to_dense_batch(graph_reps["re_mid"][1], batch.batch)
        batch_x_2, batch_mask_2 = to_dense_batch(graph_reps["re_mid"][2], batch.batch)

        graph_mids_0=self.aggrs[1](batch_x_0, batch_mask_0)
        graph_mids_1 = self.aggrs[1](batch_x_1, batch_mask_1)
        graph_mids_2 = self.aggrs[1](batch_x_2, batch_mask_2)

        out_mid = [graph_mids_0[0],graph_mids_1[0],graph_mids_2[0]]
        g_mid, w_mid = self.selector_mid(out_mid)
        return g_mid, w_mid, graph_mids_0, graph_mids_1, graph_mids_2
    
    def get_graphs_high(self,graph_reps,batch):
        batch_x_0, batch_mask_0 = to_dense_batch(graph_reps["re_high"][0], batch.batch)
        batch_x_1, batch_mask_1 = to_dense_batch(graph_reps["re_high"][1], batch.batch)
        batch_x_2, batch_mask_2 = to_dense_batch(graph_reps["re_high"][2], batch.batch)

        graph_highs_0=self.aggrs[2](batch_x_0, batch_mask_0)
        graph_highs_1 = self.aggrs[2](batch_x_1, batch_mask_1)
        graph_highs_2 = self.aggrs[2](batch_x_2, batch_mask_2)

        out_high = [graph_highs_0[0],graph_highs_1[0],graph_highs_2[0]]
        g_high, w_high = self.selector_high(out_high)
        return g_high, w_high, graph_highs_0, graph_highs_1, graph_highs_2

    
    def forward(self, data, channel_idx=-1):
        graph_reps = self.gnn(data.x, data.pe, data.edge_index, data.edge_attr, data.batch)
        g_low, w_low, graph_lows_0,graph_lows_1 = self.get_graphs_low(graph_reps,data)
        g_mid, w_mid, graph_mids_0, graph_mids_1,graph_mids_2= self.get_graphs_mid(graph_reps,data)
        g_high, w_high, graph_highs_0, graph_highs_1, graph_highs_2 = self.get_graphs_high(graph_reps,data)
        g_stack = torch.stack([g_low, g_mid, g_high], dim=1)
        score = self.attn(g_stack)
        attn_weight = torch.softmax(score, dim=1)
        g_fused = torch.sum(attn_weight * g_stack, dim=1)
        output = {}
        output['predict'] = self.pred_head(g_fused)
        return output
        
    
