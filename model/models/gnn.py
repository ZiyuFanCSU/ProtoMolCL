from typing import Any, Dict, Optional
from collections import defaultdict
import math
import copy
import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch_geometric.utils import add_self_loops, degree, softmax
from torch_geometric.nn import global_add_pool, global_mean_pool, global_max_pool
from torch_geometric.nn import GINEConv, GPSConv, MessagePassing
from torch_geometric.nn.attention import PerformerAttention
from torch_geometric.utils import to_dense_batch
from torch_geometric.nn import global_add_pool, global_mean_pool, GraphMultisetTransformer
from model.models.aggr import PromptAggr
num_atom_type = 120  
num_chirality_tag = 3
num_bond_type = 6  
num_bond_direction = 3

class GPS(torch.nn.Module):
    def __init__(self, channels: int, pe_dim: int, node_dim: int, edge_dim: int, heads: int,
                 num_layers: int, attn_type: str, dropout: float, attn_dropout: float):
        super().__init__()

        self.name = 'gps'
        self.node_emb = nn.Linear(node_dim, channels - pe_dim)
        self.pe_lin = nn.Linear(20, pe_dim)
        self.pe_norm = nn.BatchNorm1d(20)
        self.edge_emb = nn.Linear(edge_dim, channels)
        self.dropout = dropout
        self.attn_dropout = attn_dropout

        self.convs = nn.ModuleList()
        for _ in range(num_layers):
            layer = nn.Sequential(
                nn.Linear(channels, channels),
                nn.ReLU(),
                nn.Linear(channels, channels),
            )
            conv = GPSConv(channels, GINEConv(layer), heads=heads, dropout=self.dropout,
                           attn_type=attn_type, attn_kwargs={'dropout': self.attn_dropout})
            self.convs.append(conv)

    def forward(self, x, pe, edge_index, edge_attr, batch):
        x, edge_attr = x.float(), edge_attr.float()
        x_pe = self.pe_norm(pe)
        x = torch.cat((self.node_emb(x), self.pe_lin(x_pe)), 1)
        edge_attr = self.edge_emb(edge_attr)

        for i, conv in enumerate(self.convs):
            x = conv(x, edge_index, batch, edge_attr=edge_attr)

        out = global_add_pool(x, batch)
        return out, x


class GPS_multilayer(torch.nn.Module):
    def __init__(self, channels: int, pe_dim: int, node_dim: int, edge_dim: int, heads: int,
                 num_layers_low: int,num_layers_mid: int,num_layers_high: int, 
                 attn_type: str, dropout: float, attn_dropout: float):
        super().__init__()

        self.name = 'gps_multi'
        self.node_emb = nn.Linear(node_dim, channels - pe_dim)
        self.pe_lin = nn.Linear(20, pe_dim)
        self.pe_norm = nn.BatchNorm1d(20)
        self.edge_emb = nn.Linear(edge_dim, channels)
        self.dropout = dropout
        self.attn_dropout = attn_dropout

        self.convs_low = nn.ModuleList()
        for _ in range(num_layers_low):
            layer = nn.Sequential(
                nn.Linear(channels, channels),
                nn.ReLU(),
                nn.Linear(channels, channels),
            )
            conv = GPSConv(channels, GINEConv(layer), heads=heads, dropout=self.dropout,
                           attn_type=attn_type, attn_kwargs={'dropout': self.attn_dropout})
            self.convs_low.append(conv)

        self.convs_mid = nn.ModuleList()
        for _ in range(num_layers_mid):
            layer = nn.Sequential(
                nn.Linear(channels, channels),
                nn.ReLU(),
                nn.Linear(channels, channels),
            )
            conv = GPSConv(channels, GINEConv(layer), heads=heads, dropout=self.dropout,
                           attn_type=attn_type, attn_kwargs={'dropout': self.attn_dropout})
            self.convs_mid.append(conv)

        self.convs_high = nn.ModuleList()
        for _ in range(num_layers_high):
            layer = nn.Sequential(
                nn.Linear(channels, channels),
                nn.ReLU(),
                nn.Linear(channels, channels),
            )
            conv = GPSConv(channels, GINEConv(layer), heads=heads, dropout=self.dropout,
                           attn_type=attn_type, attn_kwargs={'dropout': self.attn_dropout})
            self.convs_high.append(conv)
        
        self.prompt_token = ['low','mid','high']
        self.aggrs = nn.ModuleList([PromptAggr(emb_dim=channels,
                                               num_heads=3,
                                               dropout=0,
                                               local_loss=i==(len(self.prompt_token)-1),
                                               layer_norm_out=True)
                                    for i in range(len(self.prompt_token))])

    def forward(self, x_, pe, edge_index, edge_attr, batch):
        x0, edge_attr = x_.float(), edge_attr.float()
        x_pe = self.pe_norm(pe)
        x0 = torch.cat((self.node_emb(x0), self.pe_lin(x_pe)), 1)
        edge_attr = self.edge_emb(edge_attr)

        re_low = []
        x=x0
        for i, conv in enumerate(self.convs_low):
            x = conv(x, edge_index, batch, edge_attr=edge_attr)
            re_low.append(x)

        re_mid = []
        x=x0
        for i, conv in enumerate(self.convs_mid):
            x = conv(x, edge_index, batch, edge_attr=edge_attr)
            if i>=1:
                # out_mid.append(global_add_pool(x, batch))
                re_mid.append(x)

        re_high = []
        x=x0
        for i, conv in enumerate(self.convs_high):
            x = conv(x, edge_index, batch, edge_attr=edge_attr)
            if i>=2:
                # out_high.append(global_add_pool(x, batch))
                re_high.append(x)

        return {
            "re_low": re_low,
            "re_mid": re_mid,
            "re_high": re_high,

        }