import copy
import json
import os
import torch
import pickle
import collections
import math
import lmdb
import pandas as pd
import numpy as np
from tqdm import tqdm
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs
from torch.utils import data
from torch.utils.data import Dataset
from torch_geometric.data import Data, InMemoryDataset, Batch
import torch_geometric.transforms as T
from model.utils.data import *
from model.utils.descriptors.rdNormalizedDescriptors import RDKit2DNormalized

generator = RDKit2DNormalized()

class MoleculeDataset(Dataset):
    def __init__(self, data_dir, num_candidates=5, feat_type='super_rich'):
        self.data_dir = data_dir
        self.num_candidates = num_candidates
        self.feat_type = feat_type
        self.transform = T.AddRandomWalkPE(walk_length=20, attr_name='pe')
        with open(os.path.join('./data/pretrain', 'atom_vocab_pretrain.pkl'), 'rb') as f:
            self.atom_vocab_stoi, self.atom_vocab_itos = pickle.load(f)
        with open(os.path.join("./processed_dataset/", 'precompute/data_shard_ids.pkl'), 'rb') as f:
            self.data_shard_ids = pickle.load(f)
        self.keys = list(self.data_shard_ids.keys())
        
    def get_processed_data(self, key):
        lmdb_path = os.path.join(self.data_dir, 'precompute/{}.lmdb'.format(self.data_shard_ids[key]))
        read_env = lmdb.open(lmdb_path, readonly=True, lock=False, readahead=False, meminit=False)
        with read_env.begin(write=False) as txn:
            processed_data = txn.get(key.encode())
            processed_data = json.loads(processed_data.decode())
        return processed_data

    def subgraph_mask(self, data, mol, percent=0.15, extract_context=False):
        node_num = data.x.size(0)
        label = [0] * mol.GetNumAtoms()
        mask_inds = set()
        mask_edge_inds = []
        edge_index_t = data.edge_index.transpose(0, 1)

        while len(mask_inds) < int(node_num * percent):
            center_idx = np.random.choice(list(set(np.arange(node_num)) - set(mask_inds)))
            remove_inds = [center_idx]
            while remove_inds and len(mask_inds) < int(node_num * percent):
                idx = int(remove_inds.pop(0))
                center_atom = mol.GetAtomWithIdx(idx)
                for n_atom in center_atom.GetNeighbors():
                    n_idx = n_atom.GetIdx()
                    if n_idx not in mask_inds:
                        remove_inds.append(n_idx)
                        mask_inds.add(n_idx)

                    e_idx1 = ((edge_index_t == torch.Tensor([center_idx, n_idx])).sum(dim=-1) == 2).nonzero()
                    e_idx2 = ((edge_index_t == torch.Tensor([n_idx, center_idx])).sum(dim=-1) == 2).nonzero()
                    if len(e_idx1) and len(e_idx2):
                        mask_edge_inds += [e_idx1.item(), e_idx2.item()]

                if extract_context: 
                    label[center_idx] = self.atom_vocab_stoi.get(atom_to_vocab(mol, center_atom),
                                                                 self.atom_vocab_stoi['<other>'])
                    break

            mask_inds.add(center_idx)

        for idx in mask_inds:
            if self.feat_type == 'basic':
                mask_token = torch.tensor([119, 0], dtype=torch.long)
            elif self.feat_type == 'rich':
                mask_token = (
                    onek_encoding_unk(119, list(range(1, 101))) +
                    onek_encoding_unk(0, list(range(11))) +
                    onek_encoding_unk(0, list(range(-2, 3))) +
                    onek_encoding_unk(0, list(range(5))) +
                    onek_encoding_unk(
                        Chem.rdchem.HybridizationType.UNSPECIFIED,
                        [
                            Chem.rdchem.HybridizationType.SP,
                            Chem.rdchem.HybridizationType.SP2,
                            Chem.rdchem.HybridizationType.SP3,
                            Chem.rdchem.HybridizationType.SP3D,
                            Chem.rdchem.HybridizationType.SP3D2
                        ]
                    ) +
                    [0] +
                    onek_encoding_unk(0, list(range(5))) +
                    [0] +
                    [0, 0, 0] +
                    [0]
                )

                mask_token = torch.tensor(mask_token, dtype=torch.long)
                assert mask_token.numel() == 143, \
                    f"Rich mask token dimension error: {mask_token.numel()}"
            elif self.feat_type == 'super_rich':
                mask_token = onek_encoding_unk(119, list(range(119))) + \
                             onek_encoding_unk(0, [0, 1, 2, 3, 4, 5]) + \
                             onek_encoding_unk(0, [-1, -2, 1, 2, 0]) + \
                             onek_encoding_unk(0, [0, 1, 2, 3]) + \
                             onek_encoding_unk(0, [0, 1, 2, 3, 4]) + \
                             onek_encoding_unk(0, [0, 1, 2, 3, 4]) + \
                             [0, 0] + \
                             onek_encoding_unk(0, [0, 1, 2, 3, 4, 5, 6]) + \
                             [False] * 10
                mask_token = torch.tensor(mask_token, dtype=torch.long)
            else:
                raise NotImplementedError

            data.x[idx] = mask_token

        # Apply edge attr masking
        if self.feat_type == 'basic':
            mask_token = torch.tensor([5, 0], dtype=torch.long)
        elif self.feat_type == 'super_rich':
            mask_token = onek_encoding_unk(None, [Chem.rdchem.BondType.SINGLE,
                                                  Chem.rdchem.BondType.DOUBLE,
                                                  Chem.rdchem.BondType.TRIPLE,
                                                  Chem.rdchem.BondType.AROMATIC]) + \
                         onek_encoding_unk(0, list(range(6))) + [0, 0]
            mask_token = torch.tensor(mask_token, dtype=torch.long)
        elif self.feat_type == 'rich':
            mask_token = (
                onek_encoding_unk(
                    None,
                    [
                        Chem.rdchem.BondType.SINGLE,
                        Chem.rdchem.BondType.DOUBLE,
                        Chem.rdchem.BondType.TRIPLE,
                        Chem.rdchem.BondType.AROMATIC
                    ]
                )
                +
                # Is conjugated
                [0]
                +
                # Is in ring
                [0]
                +
                # Bond stereo：None编码到unknown位，共7维
                onek_encoding_unk(
                    None,
                    [
                        Chem.rdchem.BondStereo.STEREONONE,
                        Chem.rdchem.BondStereo.STEREOANY,
                        Chem.rdchem.BondStereo.STEREOZ,
                        Chem.rdchem.BondStereo.STEREOE,
                        Chem.rdchem.BondStereo.STEREOCIS,
                        Chem.rdchem.BondStereo.STEREOTRANS
                    ]
                )
            )

            mask_token = torch.tensor(mask_token, dtype=torch.long)

            assert mask_token.numel() == 14, \
                f"Rich edge mask token dimension error: {mask_token.numel()}"
        else:
            raise NotImplementedError

        data.edge_attr[mask_edge_inds] = mask_token

        if extract_context:
            data.context_label = label
        return mask_inds

    
    def get_scaffold_index(self, mol, scaff_smi):
        index = mol.GetSubstructMatch(Chem.MolFromSmarts(scaff_smi))
        if not index:
            return list(np.arange(len(mol.GetAtoms())))
        return index

    
    def get_graph_data(self, mol):
        data = eval('mol_to_graph_data_obj_{}'.format(self.feat_type))(mol)
        data = self.transform(data)
        return data

    
    def initialize(self, data, smi, mol_fps=[], mol_mds=[], layer1=[], layer3=[], layer5=[],context_label=[]):
        data.smi = smi
        data.mol_fps = mol_fps
        data.mol_mds = mol_mds
        data.layer1 = layer1
        data.layer3 = layer3
        data.layer5 = layer5
        data.context_label = context_label
    
    
    def __len__(self):
        return len(self.keys)

    
    def __getitem__(self, idx):
        key = self.keys[idx]

        pdata = self.get_processed_data(key)
                
        smi = pdata['smi']
        ori_mol = Chem.MolFromSmiles(smi)
        pdata_md = np.array(list(generator.process(smi)), dtype=float)[1:] # 那些分子官能团和属性
        data = self.get_graph_data(ori_mol)
        mol_fp = get_fp(ori_mol)

        # sample perturbations
        si_perturb_samples = []  # scaffold-invariant perturbation samples
        num_mutations = self.num_candidates - 1

        if len(pdata['mutate']) >= num_mutations:
            si_perturb_samples += list(np.random.choice(pdata['mutate'], size=num_mutations, replace=False))
        else:
            si_perturb_samples += pdata['mutate']
            si_perturb_samples += [smi] * (num_mutations-len(pdata['mutate']))
            np.random.shuffle(si_perturb_samples)

        # get graph objects:
        candidates = []
        # original smiles
        self.initialize(data,
                        smi=smi,
                        mol_fps=[mol_fp],
                        mol_mds=[pdata_md],
                        layer1 = pdata["layer1"],
                        layer3 = pdata["layer3"],
                        layer5 = pdata["layer5"]
                        )
        
        candidates += [data]

        # random-perturbation smiles / random subgraph masking
        for i in range(5):
            data_mask = data.clone()
            self.initialize(data_mask, smi=smi)
            self.subgraph_mask(data_mask, ori_mol, extract_context=False, percent=0.15)
            candidates.append(data_mask)

        # scaffold-invariant-perturbation smiles
        for i in range(len(si_perturb_samples)):
            mol = Chem.MolFromSmiles(si_perturb_samples[i])
            data_perturb = self.get_graph_data(mol)
            mol_fp2 = get_fp(mol)
            self.initialize(data_perturb,
                            smi=si_perturb_samples[i],
                            mol_fps=[mol_fp2],
                            )                            
            candidates.append(data_perturb)

        data_context = data.clone()
        self.initialize(data_context, smi=smi,mol_fps=[mol_fp])
        self.subgraph_mask(data_context, ori_mol, extract_context=True, percent=0.15)
        candidates.append(data_context)
        
        return candidates
    
