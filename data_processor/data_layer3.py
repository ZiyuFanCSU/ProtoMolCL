import numpy as np
import pandas as pd
import json

npy_file = "./layer3_proto_k_results/K_108/molecule_soft_prototype_labels_FULL.npy"
output_csv = "./layer3.csv"
probs = np.load(npy_file)


def balance_topk_probability(
    probs,
    topk=10,
    eps=1e-12,
    T_smooth=10.0,
    T_sharp=0.5,
    entropy_low=0.45,
    entropy_high=0.85,
    alpha_smooth=0.40,
    alpha_middle=0.10,
    alpha_sharp=0.0
):
    probs = np.asarray(probs, dtype=np.float64)
    N, C = probs.shape

    if topk > C:
        raise ValueError(f"topk={topk} > C={C}")

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
        entropy = -np.sum(p_top * np.log(p_top + eps)) / np.log(topk)

        if entropy < entropy_low:
            T = T_smooth
            alpha = alpha_smooth
        elif entropy > entropy_high:
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
        entropy_list.append(entropy)
        temp_list.append(T)
        alpha_list.append(alpha)

    return new_probs, np.array(entropy_list), np.array(temp_list), np.array(alpha_list)


def print_prob_comparison(probs_before, probs_after, df, num_samples=10, topk=10):

    for i in range(min(num_samples, probs_before.shape[0])):
        p_before = probs_before[i]
        p_after = probs_after[i]

        before_top_idx = np.argsort(p_before)[-topk:][::-1]
        after_top_idx = np.argsort(p_after)[-topk:][::-1]


def print_global_stats(probs_before, probs_after):

    before_max = probs_before.max(axis=1)
    after_max = probs_after.max(axis=1)

    before_nonzero = np.sum(probs_before > 0, axis=1)
    after_nonzero = np.sum(probs_after > 0, axis=1)


new_probs, entropy_list, temp_list, alpha_list = balance_topk_probability(
    probs,
    topk=10,
    T_smooth=10.0,
    T_sharp=0.35,
    entropy_low=0.45,
    entropy_high=0.55,
    alpha_smooth=0.40,
    alpha_middle=0.10,
    alpha_sharp=0.0
)

print_prob_comparison(
    probs_before=probs,
    probs_after=new_probs,
    df=df,
    num_samples=10,
    topk=10
)

print_global_stats(
    probs_before=probs,
    probs_after=new_probs
)