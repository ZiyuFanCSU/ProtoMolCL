import numpy as np
import pandas as pd
import json


npy_file = "./layer1_proto_k_results/final_proto_soft_labels.npy"
output_csv = "./layer1.csv"
log_file = "./layer1.txt"


probs = np.load(npy_file)

def balance_topk_probability_log(probs,
                                  topk=5,
                                  T_smooth=10.0,
                                  T_sharp=0.5,
                                  entropy_low=0.45,
                                  entropy_high=0.85,
                                  alpha_smooth=0.40,
                                  alpha_middle=0.10,
                                  alpha_sharp=0.0,
                                  log_path=None):
    probs = np.asarray(probs, dtype=np.float64)
    N, C = probs.shape
    new_probs = np.zeros_like(probs)

    log_lines = []

    for i in range(N):
        p = probs[i]
        top_idx = np.argsort(p)[-topk:]
        p_top = p[top_idx].astype(np.float64)

        if p_top.sum() <= 1e-12:
            p_new_top = np.ones(topk) / topk
            new_probs[i, top_idx] = p_new_top
            continue

        p_top = p_top / (p_top.sum() + 1e-12)
        ent = -np.sum(p_top * np.log(p_top + 1e-12)) / np.log(topk)

        if ent < entropy_low:
            T = T_smooth
            alpha = alpha_smooth
        elif ent > entropy_high:
            T = T_sharp
            alpha = alpha_sharp
        else:
            T = 1.5
            alpha = alpha_middle

        p_temp = np.power(p_top + 1e-12, 1.0 / T)
        p_temp = p_temp / (p_temp.sum() + 1e-12)

        uniform_top = np.ones(topk) / topk
        p_new_top = (1.0 - alpha) * p_temp + alpha * uniform_top
        p_new_top = p_new_top / (p_new_top.sum() + 1e-12)

        new_probs[i, top_idx] = p_new_top

        if log_path:
            log_lines.append(f"Index {i} | top-{topk} positions {top_idx.tolist()}\n")
            log_lines.append(f"Before: {p[top_idx].tolist()}\n")
            log_lines.append(f"After : {p_new_top.tolist()}\n\n")

    if log_path:
        with open(log_path, "w") as f:
            f.writelines(log_lines)

    return new_probs


new_probs = balance_topk_probability_log(
    probs,
    topk=5,
    T_smooth=10.0,
    T_sharp=0.5,
    entropy_low=0.45,
    entropy_high=0.85,
    alpha_smooth=0.30,
    alpha_middle=0.10,
    alpha_sharp=0.0,
    log_path=log_file
)

