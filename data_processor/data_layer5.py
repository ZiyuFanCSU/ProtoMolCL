import numpy as np
import pandas as pd
import json

npy_file = "./layer5_proto_k_results/final_property_structural_proto_soft_labels.npy"
output_csv = "./layer5.csv"

probs = np.load(npy_file)

def normalized_entropy(p, eps=1e-12):
    p = np.asarray(p, dtype=np.float64)
    p = p / (p.sum() + eps)
    k = np.sum(p > 0)
    if k <= 1:
        return 0.0
    return -np.sum(p * np.log(p + eps)) / np.log(len(p))


def balance_topk_probability(
    probs,
    topk=5,
    eps=1e-12,
    T_smooth=8.0,
    T_sharp=0.5,
    entropy_low=0.45,
    entropy_high=0.85,
    alpha_smooth=0.35,
    alpha_middle=0.15,
    alpha_sharp=0.0
):

    probs = np.asarray(probs, dtype=np.float64)
    N, C = probs.shape

    new_probs = np.zeros_like(probs)
    entropy_list = []
    temp_list = []
    alpha_list = []

    for i in range(N):
        p = probs[i]
        top_idx = np.argsort(p)[-topk:]
        p_top = p[top_idx].astype(np.float64)
        if p_top.sum() <= eps:
            p_new_top = np.ones(topk) / topk
            new_probs[i, top_idx] = p_new_top
            entropy_list.append(1.0)
            temp_list.append(1.0)
            alpha_list.append(1.0)
            continue
        p_top = p_top / (p_top.sum() + eps)
        ent = -np.sum(p_top * np.log(p_top + eps)) / np.log(topk)
        if ent < entropy_low:
            T = T_smooth
            alpha = alpha_smooth
        elif ent > entropy_high:
            T = T_sharp
            alpha = alpha_sharp
        else:
            T = 1.5
            alpha = alpha_middle

        p_temp = np.power(p_top + eps, 1.0 / T)
        p_temp = p_temp / (p_temp.sum() + eps)

        uniform_top = np.ones(topk) / topk
        p_new_top = (1.0 - alpha) * p_temp + alpha * uniform_top
        p_new_top = p_new_top / (p_new_top.sum() + eps)
        new_probs[i, top_idx] = p_new_top

        entropy_list.append(ent)
        temp_list.append(T)
        alpha_list.append(alpha)

    return new_probs, np.array(entropy_list), np.array(temp_list), np.array(alpha_list)


new_probs, entropy_list, temp_list, alpha_list = balance_topk_probability(
    probs,
    topk=5,
    T_smooth=10.0,
    T_sharp=0.5,
    entropy_low=0.45,
    entropy_high=0.85,
    alpha_smooth=0.40,
    alpha_middle=0.10,
    alpha_sharp=0.0
)
