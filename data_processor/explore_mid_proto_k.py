import os
import json
import time
import argparse
import warnings
from collections import Counter
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from tqdm import tqdm
from scipy import sparse
from sklearn.cluster import MiniBatchKMeans
from sklearn.decomposition import TruncatedSVD
from sklearn.metrics import silhouette_score, calinski_harabasz_score, davies_bouldin_score
from sklearn.preprocessing import normalize
from rdkit import Chem, DataStructs
from rdkit.Chem import BRICS, AllChem
from rdkit import RDLogger

RDLogger.DisableLog("rdApp.*")
warnings.filterwarnings("ignore")

def parse_args():
    parser = argparse.ArgumentParser(
        description="BRICS-informed molecular prototype construction with similarity-smoothed fragment composition."
    )
    parser.add_argument("--csv_path", type=str, default="./zinc_200w.csv")
    parser.add_argument("--smiles_col", type=str, default="smiles")
    parser.add_argument("--out_dir", type=str, default="./layer3_proto_k_results")
    parser.add_argument("--k_min", type=int, default=100)
    parser.add_argument("--k_max", type=int, default=110)
    parser.add_argument("--batch_size", type=int, default=8192)
    parser.add_argument("--max_iter", type=int, default=100)
    parser.add_argument("--n_init", type=int, default=3)
    parser.add_argument("--radius", type=int, default=2)
    parser.add_argument("--n_bits", type=int, default=2048)
    parser.add_argument("--min_fragment_freq", type=int, default=10,
                        help="Fragments appearing fewer times than this will be ignored to reduce noise and memory.")
    parser.add_argument("--fallback_full_mol", default=True,
                        help="Use full molecule canonical SMILES as fallback when BRICS fails or gives no fragment.")
    parser.add_argument("--sim_sample_pairs", type=int, default=200000,
                        help="Number of random fragment pairs for estimating similarity distribution.")
    parser.add_argument("--sim_threshold", type=float, default=-1.0,
                        help="If >0, use this threshold directly. If <=0, auto-select threshold from sampled similarities.")
    parser.add_argument("--sim_quantile", type=float, default=0.85,
                        help="Quantile of positive sampled similarities for auto threshold selection.")
    parser.add_argument("--min_auto_threshold", type=float, default=0.28,
                        help="Lower bound for auto-selected similarity threshold.")
    parser.add_argument("--max_auto_threshold", type=float, default=0.75,
                        help="Upper bound for auto-selected similarity threshold.")
    parser.add_argument("--max_sim_neighbors", type=int, default=200,
                        help="For each fragment, keep at most this many neighbors above threshold.")
    parser.add_argument("--sim_chunk_size", type=int, default=256)
    parser.add_argument("--use_tfidf", default=True,
                        help="Use TF-IDF weighting before similarity smoothing.")
    parser.add_argument("--svd_dim", type=int, default=256,
                        help="Reduce smoothed sparse vectors by TruncatedSVD before clustering. Set <=0 to disable.")
    parser.add_argument("--svd_sample_size", type=int, default=300000,
                        help="Fit SVD on a sample of molecules to save memory/time.")
    parser.add_argument("--save_soft_labels", default=True,
                        help="Save molecule-level soft prototype affiliation vectors for each K.")
    parser.add_argument("--soft_tau", type=float, default=1.0,
                        help="Temperature for converting distances to soft prototype probabilities.")
    parser.add_argument("--soft_batch_size", type=int, default=10000,
                        help="Batch size for computing soft prototype probabilities.")
    parser.add_argument("--metric_sample_size", type=int, default=50000,
                        help="Sample size for clustering validity metrics.")
    parser.add_argument("--top_fragments", type=int, default=20)
    parser.add_argument("--seed", type=int, default=42)

    return parser.parse_args()


def ensure_dir(path):
    os.makedirs(path, exist_ok=True)


def canonical_smiles(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, None
    return Chem.MolToSmiles(mol, canonical=True), mol


def canonical_fragment_smi(frag_smi):
    mol = Chem.MolFromSmiles(frag_smi)
    if mol is None:
        return None
    return Chem.MolToSmiles(mol, canonical=True)


def brics_fragments_from_mol(mol, fallback_smi=None):
    try:
        broken = BRICS.BreakBRICSBonds(mol)
        frag_mols = Chem.GetMolFrags(broken, asMols=True, sanitizeFrags=True)
        frags = []
        for fm in frag_mols:
            fs = Chem.MolToSmiles(fm, canonical=True)
            cfs = canonical_fragment_smi(fs)
            if cfs is not None:
                frags.append(cfs)
        if len(frags) == 0 and fallback_smi is not None:
            frags = [fallback_smi]
        return frags
    except Exception:
        if fallback_smi is not None:
            return [fallback_smi]
        return []


def build_fragment_count_matrix(smiles_list, min_fragment_freq=10, fallback_full_mol=True):
    all_mol_frag_lists = []
    frag_global_counter = Counter()
    valid_indices = []
    valid_smiles = []

    print("[1/8] Parsing molecules and extracting BRICS fragments...")
    for idx, smi in enumerate(tqdm(smiles_list, total=len(smiles_list))):
        if not isinstance(smi, str) or len(smi.strip()) == 0:
            continue
        can_smi, mol = canonical_smiles(smi)
        if mol is None:
            continue
        fallback = can_smi if fallback_full_mol else None
        frags = brics_fragments_from_mol(mol, fallback_smi=fallback)
        if len(frags) == 0:
            continue
        all_mol_frag_lists.append(frags)
        frag_global_counter.update(frags)
        valid_indices.append(idx)
        valid_smiles.append(can_smi)

    vocab = sorted([f for f, cnt in frag_global_counter.items() if cnt >= min_fragment_freq])
    frag2id = {f: i for i, f in enumerate(vocab)}

    print(f"Valid molecules: {len(valid_smiles)}")
    print(f"Unique BRICS fragments before filtering: {len(frag_global_counter)}")
    print(f"Unique BRICS fragments after min_freq >= {min_fragment_freq}: {len(vocab)}")

    rows, cols, data = [], [], []
    print("[2/8] Building molecule-fragment count matrix...")
    for i, frags in enumerate(tqdm(all_mol_frag_lists, total=len(all_mol_frag_lists))):
        cnt = Counter([f for f in frags if f in frag2id])
        for f, c in cnt.items():
            rows.append(i)
            cols.append(frag2id[f])
            data.append(float(c))

    C = sparse.csr_matrix((data, (rows, cols)), shape=(len(valid_smiles), len(vocab)), dtype=np.float32)
    nonzero_mol_mask = np.asarray(C.sum(axis=1)).ravel() > 0
    C = C[nonzero_mol_mask]
    valid_smiles = [s for s, keep in zip(valid_smiles, nonzero_mol_mask) if keep]
    valid_indices = [i for i, keep in zip(valid_indices, nonzero_mol_mask) if keep]

    print(f"Molecules kept after fragment filtering: {C.shape[0]}")
    return C, vocab, frag_global_counter, valid_smiles, valid_indices


def fragment_fingerprints(vocab, radius=2, n_bits=2048):
    print("[3/8] Computing Morgan fingerprints for unique BRICS fragments...")
    fps = []
    kept_vocab = []
    old_to_new = {}
    for old_id, fs in enumerate(tqdm(vocab)):
        mol = Chem.MolFromSmiles(fs)
        if mol is None:
            continue
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=n_bits)
        old_to_new[old_id] = len(fps)
        fps.append(fp)
        kept_vocab.append(fs)
    if len(kept_vocab) < len(vocab):
        print(f"Warning: {len(vocab) - len(kept_vocab)} fragments failed fingerprint generation.")
    return fps, kept_vocab, old_to_new


def remap_count_matrix(C, old_to_new):
    old_cols = np.array(sorted(old_to_new.keys()), dtype=np.int64)
    return C[:, old_cols].tocsr()


def sample_similarity_distribution(fps, n_pairs=200000, seed=42):
    print("[4/8] Sampling pairwise Tanimoto similarities for threshold analysis...")
    rng = np.random.default_rng(seed)
    m = len(fps)
    if m < 2:
        raise ValueError("Need at least two fragments to compute similarities.")
    n_pairs = min(n_pairs, m * max(1, min(m - 1, 100)))
    sims = np.zeros(n_pairs, dtype=np.float32)
    for t in tqdm(range(n_pairs)):
        i = int(rng.integers(0, m))
        j = int(rng.integers(0, m - 1))
        if j >= i:
            j += 1
        sims[t] = DataStructs.TanimotoSimilarity(fps[i], fps[j])

    positive = sims[sims > 0]
    report = {
        "sample_pairs": int(n_pairs),
        "all_min": float(np.min(sims)),
        "all_mean": float(np.mean(sims)),
        "all_median": float(np.median(sims)),
        "all_p75": float(np.quantile(sims, 0.75)),
        "all_p90": float(np.quantile(sims, 0.90)),
        "all_p95": float(np.quantile(sims, 0.95)),
        "all_p99": float(np.quantile(sims, 0.99)),
        "all_max": float(np.max(sims)),
        "positive_ratio": float(len(positive) / len(sims)),
    }
    if len(positive) > 0:
        report.update({
            "pos_min": float(np.min(positive)),
            "pos_mean": float(np.mean(positive)),
            "pos_median": float(np.median(positive)),
            "pos_p75": float(np.quantile(positive, 0.75)),
            "pos_p90": float(np.quantile(positive, 0.90)),
            "pos_p95": float(np.quantile(positive, 0.95)),
            "pos_p99": float(np.quantile(positive, 0.99)),
            "pos_max": float(np.max(positive)),
        })
    return sims, report


def plot_similarity_distribution(sims, threshold, out_dir):
    os.makedirs(out_dir, exist_ok=True)

    sims = np.asarray(sims, dtype=np.float32)

    plt.figure(figsize=(7, 5))
    plt.hist(sims, bins=80, alpha=0.85)
    plt.axvline(threshold, linestyle="--", linewidth=2, label=f"threshold={threshold:.3f}")
    plt.xlabel("Pairwise Tanimoto Similarity")
    plt.ylabel("Frequency")
    plt.title("Sampled Pairwise BRICS Fragment Similarity Distribution")
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, "pairwise_tanimoto_distribution.png"), dpi=300)
    plt.close()

    positive = sims[sims > 0]
    if len(positive) > 0:
        plt.figure(figsize=(7, 5))
        plt.hist(positive, bins=80, alpha=0.85)
        plt.axvline(threshold, linestyle="--", linewidth=2, label=f"threshold={threshold:.3f}")
        plt.xlabel("Positive Pairwise Tanimoto Similarity")
        plt.ylabel("Frequency")
        plt.title("Positive Pairwise BRICS Fragment Similarity Distribution")
        plt.legend()
        plt.tight_layout()
        plt.savefig(os.path.join(out_dir, "positive_pairwise_tanimoto_distribution.png"), dpi=300)
        plt.close()

    # Also save a log-scaled version because many random fragment pairs may have zero or very low similarity.
    plt.figure(figsize=(7, 5))
    plt.hist(sims, bins=80, alpha=0.85)
    plt.axvline(threshold, linestyle="--", linewidth=2, label=f"threshold={threshold:.3f}")
    plt.yscale("log")
    plt.xlabel("Pairwise Tanimoto Similarity")
    plt.ylabel("Frequency (log scale)")
    plt.title("Sampled Pairwise Similarity Distribution (Log Scale)")
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, "pairwise_tanimoto_distribution_log.png"), dpi=300)
    plt.close()


def auto_select_threshold(report, sim_quantile=0.90, min_thr=0.35, max_thr=0.75):
    key = f"pos_p{int(sim_quantile * 100)}"
    if key in report:
        thr = report[key]
    else:
        thr = report.get("all_p95", 0.4)
    return float(max(min_thr, min(max_thr, thr)))


def build_sparse_similarity_matrix(fps, threshold=0.4, max_neighbors=200, chunk_size=256):
    print(f"[5/8] Building sparse fragment similarity matrix with threshold={threshold:.4f}...")
    m = len(fps)
    rows, cols, vals = [], [], []

    for start in tqdm(range(0, m, chunk_size)):
        end = min(m, start + chunk_size)
        for i in range(start, end):
            sim_arr = np.asarray(DataStructs.BulkTanimotoSimilarity(fps[i], fps), dtype=np.float32)
            idx = np.where(sim_arr >= threshold)[0]

            if len(idx) > max_neighbors + 1:
                self_idx = np.array([i], dtype=np.int64)
                non_self = idx[idx != i]
                if len(non_self) > max_neighbors:
                    top_local = np.argpartition(sim_arr[non_self], -max_neighbors)[-max_neighbors:]
                    non_self = non_self[top_local]
                idx = np.unique(np.concatenate([self_idx, non_self]))

            if i not in idx:
                idx = np.append(idx, i)

            rows.extend([i] * len(idx))
            cols.extend(idx.tolist())
            vals.extend(sim_arr[idx].tolist())

    S = sparse.csr_matrix((vals, (rows, cols)), shape=(m, m), dtype=np.float32)
    S = S.maximum(S.T).tocsr()
    S = normalize(S, norm="l1", axis=1, copy=False)

    nnz_per_row = np.diff(S.indptr)
    sim_report = {
        "num_fragments": int(m),
        "threshold": float(threshold),
        "max_neighbors": int(max_neighbors),
        "nnz": int(S.nnz),
        "avg_neighbors": float(nnz_per_row.mean()),
        "median_neighbors": float(np.median(nnz_per_row)),
        "min_neighbors": int(nnz_per_row.min()),
        "max_neighbors_observed": int(nnz_per_row.max()),
        "density": float(S.nnz / (m * m)),
    }
    return S, sim_report


def apply_tfidf(C):
    n = C.shape[0]
    df = np.diff(C.tocsc().indptr).astype(np.float32)
    idf = np.log((1.0 + n) / (1.0 + df)) + 1.0
    return C.multiply(idf).tocsr(), idf


def smooth_molecular_vectors(C, S, use_tfidf=True):
    print("[6/8] Constructing similarity-smoothed molecular BRICS composition vectors...")
    if use_tfidf:
        Cw, idf = apply_tfidf(C)
    else:
        Cw, idf = C, None
    X = Cw.dot(S).tocsr()
    X = normalize(X, norm="l1", axis=1, copy=False)
    return X, idf


def fit_svd_if_needed(X, dim=256, sample_size=300000, seed=42):
    if dim is None or dim <= 0:
        return X, None
    dim = min(dim, max(2, X.shape[1] - 1))
    print(f"[7/8] Reducing dimension with TruncatedSVD to {dim} dims...")
    rng = np.random.default_rng(seed)
    n = X.shape[0]
    if n > sample_size:
        sample_idx = rng.choice(n, size=sample_size, replace=False)
        X_fit = X[sample_idx]
    else:
        X_fit = X
    svd = TruncatedSVD(n_components=dim, random_state=seed)
    svd.fit(X_fit)
    X_red = svd.transform(X).astype(np.float32)
    norm = np.linalg.norm(X_red, axis=1, keepdims=True) + 1e-12
    X_red = X_red / norm
    return X_red, svd


def evaluate_clustering(X_eval, labels):
    uniq = np.unique(labels)
    if len(uniq) < 2 or len(uniq) >= len(labels):
        return {"silhouette": np.nan, "calinski_harabasz": np.nan, "davies_bouldin": np.nan}
    return {
        "silhouette": float(silhouette_score(X_eval, labels, metric="euclidean")),
        "calinski_harabasz": float(calinski_harabasz_score(X_eval, labels)),
        "davies_bouldin": float(davies_bouldin_score(X_eval, labels)),
    }


def compute_soft_labels(X, centers, tau=1.0, batch_size=10000):
    if tau <= 0:
        raise ValueError("--soft_tau must be positive.")

    centers = centers.astype(np.float32)
    center_norm = np.sum(centers ** 2, axis=1, keepdims=True).T  # [1, K]
    all_probs = []

    for start in tqdm(range(0, X.shape[0], batch_size), desc="Computing soft prototype labels"):
        end = min(start + batch_size, X.shape[0])
        Xb = X[start:end].astype(np.float32)

        # Squared Euclidean distance:
        # ||x - c||^2 = ||x||^2 + ||c||^2 - 2 x c^T
        x_norm = np.sum(Xb ** 2, axis=1, keepdims=True)  # [B, 1]
        dist = x_norm + center_norm - 2.0 * np.dot(Xb, centers.T)
        dist = np.maximum(dist, 0.0)

        logits = -dist / tau
        logits = logits - np.max(logits, axis=1, keepdims=True)
        probs = np.exp(logits)
        probs = probs / (np.sum(probs, axis=1, keepdims=True) + 1e-12)
        all_probs.append(probs.astype(np.float32))

    return np.vstack(all_probs)


def cluster_size_stats(labels):
    counts = Counter(labels.tolist())
    arr = np.array(list(counts.values()), dtype=np.float32)
    return {
        "min_cluster_size": int(arr.min()),
        "max_cluster_size": int(arr.max()),
        "mean_cluster_size": float(arr.mean()),
        "median_cluster_size": float(np.median(arr)),
        "max_cluster_ratio": float(arr.max() / arr.sum()),
        "min_cluster_ratio": float(arr.min() / arr.sum()),
    }


def top_fragments_by_cluster(C, labels, vocab, top_n=20):
    print("Generating top-fragment interpretation...")
    C_bin = C.copy()
    C_bin.data = np.ones_like(C_bin.data)
    global_presence = np.asarray(C_bin.mean(axis=0)).ravel()
    results = {}
    for k in sorted(np.unique(labels)):
        idx = np.where(labels == k)[0]
        cluster_presence = np.asarray(C_bin[idx].mean(axis=0)).ravel()
        enrich = cluster_presence - global_presence
        top_idx = np.argsort(-enrich)[:top_n]
        rows = []
        for j in top_idx:
            rows.append({
                "fragment": vocab[j],
                "cluster_presence": float(cluster_presence[j]),
                "global_presence": float(global_presence[j]),
                "presence_enrichment": float(enrich[j]),
            })
        results[int(k)] = rows
    return results


def save_cluster_interpretation_md(path, labels, top_frag_dict):
    counts = Counter(labels.tolist())
    total = len(labels)
    with open(path, "w", encoding="utf-8") as f:
        f.write("# BRICS-informed molecular prototype interpretation\n\n")
        f.write(f"Total molecules: {total}\n\n")
        f.write(f"Number of prototypes: {len(counts)}\n\n")
        f.write("## How to read this report\n")
        f.write("- Clustering objects are molecules, not fragments.\n")
        f.write("- `Top BRICS fragments` lists fragments enriched in each molecular prototype.\n")
        f.write("- `cluster_presence` means the fraction of molecules in this prototype containing the fragment.\n")
        f.write("- `presence_enrichment = cluster_presence - global_presence`.\n\n")
        for k in sorted(counts):
            f.write(f"## Prototype {k}\n")
            f.write(f"- Count: {counts[k]}\n")
            f.write(f"- Ratio: {counts[k] / total:.4f}\n\n")
            f.write("### Top BRICS fragments\n")
            for item in top_frag_dict.get(k, []):
                f.write(
                    f"- `{item['fragment']}`: "
                    f"presence={item['cluster_presence']:.4f}, "
                    f"global={item['global_presence']:.4f}, "
                    f"enrich={item['presence_enrichment']:.4f}\n"
                )
            f.write("\n")


def normalize_model_selection_scores(metrics_df):
    df = metrics_df.copy()
    def minmax(x):
        x = np.asarray(x, dtype=np.float64)
        mn, mx = np.nanmin(x), np.nanmax(x)
        if abs(mx - mn) < 1e-12:
            return np.ones_like(x) * 0.5
        return (x - mn) / (mx - mn + 1e-12)
    sil_n = minmax(df["silhouette"].values)
    ch_n = minmax(df["calinski_harabasz"].values)
    db_n = 1.0 - minmax(df["davies_bouldin"].values)
    df["silhouette_norm"] = sil_n
    df["calinski_harabasz_norm"] = ch_n
    df["davies_bouldin_inv_norm"] = db_n
    df["selection_score"] = 0.5 * sil_n + 0.3 * ch_n + 0.2 * db_n
    return df


def main():
    args = parse_args()
    ensure_dir(args.out_dir)
    t0 = time.time()

    print("Arguments:")
    print(json.dumps(vars(args), indent=2, ensure_ascii=False))

    df = pd.read_csv(args.csv_path)
    if args.smiles_col not in df.columns:
        raise ValueError(f"Column `{args.smiles_col}` not found in CSV. Available columns: {list(df.columns)}")
    smiles_list = df[args.smiles_col].tolist()

    C, vocab, frag_counter, valid_smiles, valid_indices = build_fragment_count_matrix(
        smiles_list, min_fragment_freq=args.min_fragment_freq, fallback_full_mol=args.fallback_full_mol
    )

    fps, kept_vocab, old_to_new = fragment_fingerprints(vocab, radius=args.radius, n_bits=args.n_bits)
    C = remap_count_matrix(C, old_to_new)
    vocab = kept_vocab

    sims_sample, threshold_report = sample_similarity_distribution(fps, n_pairs=args.sim_sample_pairs, seed=args.seed)
    if args.sim_threshold > 0:
        threshold = float(args.sim_threshold)
        threshold_report["threshold_source"] = "manual"
    else:
        threshold = auto_select_threshold(
            threshold_report,
            sim_quantile=args.sim_quantile,
            min_thr=args.min_auto_threshold,
            max_thr=args.max_auto_threshold,
        )
        threshold_report["threshold_source"] = "auto"
    threshold_report["selected_threshold"] = threshold
    plot_similarity_distribution(sims_sample, threshold, args.out_dir)

    with open(os.path.join(args.out_dir, "similarity_threshold_report.json"), "w", encoding="utf-8") as f:
        json.dump(threshold_report, f, indent=2, ensure_ascii=False)
    pd.DataFrame({"sampled_tanimoto": sims_sample}).to_csv(
        os.path.join(args.out_dir, "sampled_pairwise_tanimoto.csv"), index=False
    )

    S, sim_graph_report = build_sparse_similarity_matrix(
        fps, threshold=threshold, max_neighbors=args.max_sim_neighbors, chunk_size=args.sim_chunk_size
    )
    with open(os.path.join(args.out_dir, "fragment_similarity_graph_report.json"), "w", encoding="utf-8") as f:
        json.dump(sim_graph_report, f, indent=2, ensure_ascii=False)

    X_sparse, idf = smooth_molecular_vectors(C, S, use_tfidf=args.use_tfidf)

    sparse.save_npz(os.path.join(args.out_dir, "molecule_brics_count_matrix.npz"), C)
    sparse.save_npz(os.path.join(args.out_dir, "fragment_similarity_matrix.npz"), S)
    sparse.save_npz(os.path.join(args.out_dir, "molecule_brics_smoothed_matrix.npz"), X_sparse)
    with open(os.path.join(args.out_dir, "brics_vocab.txt"), "w", encoding="utf-8") as f:
        for frag in vocab:
            f.write(frag + "\n")
    pd.DataFrame({"valid_original_index": valid_indices, "canonical_smiles": valid_smiles}).to_csv(
        os.path.join(args.out_dir, "valid_molecules.csv"), index=False
    )

    X_cluster, svd_model = fit_svd_if_needed(
        X_sparse, dim=args.svd_dim, sample_size=args.svd_sample_size, seed=args.seed
    )

    rng = np.random.default_rng(args.seed)
    n = X_cluster.shape[0]
    metric_n = min(args.metric_sample_size, n)
    metric_idx = rng.choice(n, size=metric_n, replace=False)

    metrics_rows = []
    print("[8/8] Clustering molecules for K range...")
    for K in range(args.k_min, args.k_max + 1):
        print(f"\nRunning MiniBatchKMeans with K={K}...")
        k_dir = os.path.join(args.out_dir, f"K_{K}")
        ensure_dir(k_dir)

        km = MiniBatchKMeans(
            n_clusters=K,
            random_state=args.seed,
            batch_size=args.batch_size,
            max_iter=args.max_iter,
            n_init=args.n_init,
            verbose=0,
        )
        labels = km.fit_predict(X_cluster)

        metric = evaluate_clustering(X_cluster[metric_idx], labels[metric_idx])
        size_stat = cluster_size_stats(labels)
        row = {"K": K, **metric, **size_stat}
        metrics_rows.append(row)

        pd.DataFrame({
            "valid_original_index": valid_indices,
            "canonical_smiles": valid_smiles,
            "prototype_label": labels,
        }).to_csv(os.path.join(k_dir, "molecule_prototype_labels.csv"), index=False)
        np.save(os.path.join(k_dir, "cluster_centers.npy"), km.cluster_centers_)

        if args.save_soft_labels:
            soft_probs = compute_soft_labels(
                X_cluster,
                km.cluster_centers_,
                tau=args.soft_tau,
                batch_size=args.soft_batch_size,
            )

            np.save(os.path.join(k_dir, "molecule_soft_prototype_labels.npy"), soft_probs)
            with open(os.path.join(k_dir, "soft_label_info.json"), "w", encoding="utf-8") as f:
                json.dump({
                    "soft_label_file_npy": "molecule_soft_prototype_labels.npy",
                    "soft_label_file_csv": "molecule_soft_prototype_labels.csv",
                    "shape": [int(soft_probs.shape[0]), int(soft_probs.shape[1])],
                    "tau": float(args.soft_tau),
                    "definition": "p_ik = softmax_k(-||x_i - mu_k||_2^2 / tau)",
                    "note": "Rows are aligned with valid_molecules.csv and molecule_prototype_labels.csv."
                }, f, indent=2, ensure_ascii=False)

            del soft_probs

        top_frag_dict = top_fragments_by_cluster(C, labels, vocab, top_n=args.top_fragments)
        with open(os.path.join(k_dir, "top_brics_fragments.json"), "w", encoding="utf-8") as f:
            json.dump(top_frag_dict, f, indent=2, ensure_ascii=False)
        save_cluster_interpretation_md(os.path.join(k_dir, "cluster_interpretation.md"), labels, top_frag_dict)
        with open(os.path.join(k_dir, "metrics.json"), "w", encoding="utf-8") as f:
            json.dump(row, f, indent=2, ensure_ascii=False)

        print(json.dumps(row, indent=2, ensure_ascii=False))

    metrics_df = pd.DataFrame(metrics_rows)
    metrics_df = normalize_model_selection_scores(metrics_df)
    metrics_df = metrics_df.sort_values("selection_score", ascending=False)
    metrics_df.to_csv(os.path.join(args.out_dir, "k_selection_metrics.csv"), index=False)

    best_K = int(metrics_df.iloc[0]["K"])
    summary = {
        "best_K_by_selection_score": best_K,
        "selection_metrics_csv": "k_selection_metrics.csv",
        "elapsed_seconds": float(time.time() - t0),
        "note": (
            "The best K is selected by normalized Silhouette, CH, and inverse DB scores. "
            "Please also inspect cluster_interpretation.md under each K folder to avoid "
            "fragment-count or molecular-size shortcuts."
        ),
    }
    with open(os.path.join(args.out_dir, "run_summary.json"), "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2, ensure_ascii=False)

    print("\nDone.")
    print(json.dumps(summary, indent=2, ensure_ascii=False))
    print(f"Results saved to: {args.out_dir}")


if __name__ == "__main__":
    main()
