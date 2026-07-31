import os
import argparse
import warnings
import numpy as np
import pandas as pd
from tqdm import tqdm
from sklearn.cluster import MiniBatchKMeans
from sklearn.metrics import silhouette_score, calinski_harabasz_score, davies_bouldin_score
from data_process2.descriptors.rdNormalizedDescriptors import RDKit2DNormalized


def compute_rdkit2d_vectors(csv_path, smiles_col="smiles", out_dir="./proto_k_results"):
    os.makedirs(out_dir, exist_ok=True)

    df = pd.read_csv(csv_path)

    if smiles_col not in df.columns:
        raise ValueError(f"Column not found in CSV: {smiles_col}. Available columns: {list(df.columns)}")

    generator = RDKit2DNormalized()
    smiles_list = df[smiles_col].astype(str).tolist()

    valid_smiles = []
    valid_indices = []
    features = []
    failed = []

    for idx, smi in tqdm(list(enumerate(smiles_list)), desc="Computing RDKit2DNormalized"):
        try:
            vec = np.array(list(generator.process(smi)), dtype=np.float32)
            vec = vec[1:]

            if np.any(np.isnan(vec)) or np.any(np.isinf(vec)):
                failed.append((idx, smi, "nan_or_inf"))
                continue

            valid_smiles.append(smi)
            valid_indices.append(idx)
            features.append(vec)

        except Exception as e:
            failed.append((idx, smi, str(e)))
            continue

    if len(features) == 0:
        raise RuntimeError("No molecular features were computed successfully. Please check the SMILES column and the RDKit2DNormalized environment.")

    X = np.vstack(features).astype(np.float32)

    print(f"\nNumber of valid molecules: {len(valid_smiles)}")
    print(f"Number of failed molecules: {len(failed)}")
    print(f"Feature dimension: {X.shape[1]}")
    print(f"Feature dtype: {X.dtype}")

    np.save(os.path.join(out_dir, "rdkit2d_normalized_features.npy"), X)
    valid_df = df.iloc[valid_indices].copy()
    valid_df["valid_index"] = np.arange(len(valid_df))
    valid_df.to_csv(os.path.join(out_dir, "valid_smiles.csv"), index=False)

    if len(failed) > 0:
        failed_df = pd.DataFrame(failed, columns=["row_index", "smiles", "reason"])
        failed_df.to_csv(os.path.join(out_dir, "failed_smiles.csv"), index=False)

    return X, valid_df


def sample_for_metric(X, labels, sample_size=20000, random_state=42):
    n = X.shape[0]
    if n <= sample_size:
        return X, labels

    rng = np.random.default_rng(random_state)
    idx = rng.choice(n, size=sample_size, replace=False)
    return X[idx], labels[idx]


def safe_metric_scores(X_metric, labels_metric, random_state=42):
    unique_labels = np.unique(labels_metric)
    if len(unique_labels) < 2:
        return np.nan, np.nan, np.nan

    sil = silhouette_score(X_metric, labels_metric, random_state=random_state)
    ch = calinski_harabasz_score(X_metric, labels_metric)
    db = davies_bouldin_score(X_metric, labels_metric)
    return sil, ch, db


def explore_kmeans_k(
    X,
    k_min=2,
    k_max=20,
    out_dir="./proto_k_results",
    random_state=42,
    min_cluster_ratio=0.005,
    min_cluster_size=200,
    batch_size=8192,
    metric_sample_size=20000,
):

    os.makedirs(out_dir, exist_ok=True)

    results = []
    n_samples = X.shape[0]

    for k in tqdm(range(k_min, k_max + 1), desc="Exploring K"):
        if k >= n_samples:
            continue
        try:
            kmeans = MiniBatchKMeans(
                n_clusters=k,
                random_state=random_state,
                batch_size=batch_size,
                n_init=10,
                max_iter=300,
                reassignment_ratio=0.01,
                verbose=0,
            )

            labels = kmeans.fit_predict(X)

            counts = np.bincount(labels, minlength=k)
            min_count = counts.min()

            required_min_size = max(min_cluster_size, int(n_samples * min_cluster_ratio))
            valid_cluster_size = min_count >= required_min_size

            X_metric, labels_metric = sample_for_metric(
                X,
                labels,
                sample_size=metric_sample_size,
                random_state=random_state,
            )

            sil, ch, db = safe_metric_scores(
                X_metric,
                labels_metric,
                random_state=random_state,
            )

            inertia = kmeans.inertia_

            results.append({
                "K": k,
                "silhouette_sample": sil,              
                "calinski_harabasz_sample": ch,        
                "davies_bouldin_sample": db,           
                "inertia": inertia,                    
                "min_cluster_size": int(min_count),
                "max_cluster_size": int(counts.max()),
                "required_min_size": int(required_min_size),
                "valid_cluster_size": valid_cluster_size,
                "cluster_counts": counts.tolist(),
            })

        except Exception as e:
            warnings.warn(f"K={k} clustering failed: {e}")

    result_df = pd.DataFrame(results)
    if len(result_df) == 0:
        raise RuntimeError("No clustering run was completed successfully for any K. Please check the input features.")

    result_df.to_csv(os.path.join(out_dir, "kmeans_k_search_all.csv"), index=False)
    valid_df = result_df[result_df["valid_cluster_size"] == True].copy()
    if len(valid_df) == 0:
        print("\nWarning: No K satisfies the minimum cluster-size requirement. Selection will be made from all K values.")
        valid_df = result_df.copy()

    metric_cols = ["silhouette_sample", "calinski_harabasz_sample", "davies_bouldin_sample"]
    valid_df = valid_df.dropna(subset=metric_cols).copy()

    if len(valid_df) == 0:
        raise RuntimeError("All clustering metrics are NaN for every K. Please check the clustering results.")

    eps = 1e-12

    valid_df["sil_norm"] = (
        valid_df["silhouette_sample"] - valid_df["silhouette_sample"].min()
    ) / (valid_df["silhouette_sample"].max() - valid_df["silhouette_sample"].min() + eps)

    valid_df["ch_norm"] = (
        valid_df["calinski_harabasz_sample"] - valid_df["calinski_harabasz_sample"].min()
    ) / (valid_df["calinski_harabasz_sample"].max() - valid_df["calinski_harabasz_sample"].min() + eps)

    valid_df["db_norm"] = (
        valid_df["davies_bouldin_sample"].max() - valid_df["davies_bouldin_sample"]
    ) / (valid_df["davies_bouldin_sample"].max() - valid_df["davies_bouldin_sample"].min() + eps)

    valid_df["score"] = (
        0.5 * valid_df["sil_norm"]
        + 0.3 * valid_df["ch_norm"]
        + 0.2 * valid_df["db_norm"]
    )

    best_row = valid_df.sort_values("score", ascending=False).iloc[0]
    best_k = int(best_row["K"])

    valid_df.to_csv(os.path.join(out_dir, "kmeans_k_search_valid_scored.csv"), index=False)

    print("\n========== K SEARCH RESULTS ==========")
    print(result_df[[
        "K",
        "silhouette_sample",
        "calinski_harabasz_sample",
        "davies_bouldin_sample",
        "min_cluster_size",
        "max_cluster_size",
        "required_min_size",
        "valid_cluster_size",
    ]])

    print("\n========== RECOMMENDED K ==========")
    print(f"Best K = {best_k}")
    print(f"Silhouette(sample) = {best_row['silhouette_sample']:.4f}")
    print(f"Calinski-Harabasz(sample) = {best_row['calinski_harabasz_sample']:.4f}")
    print(f"Davies-Bouldin(sample) = {best_row['davies_bouldin_sample']:.4f}")
    print(f"Min cluster size = {int(best_row['min_cluster_size'])}")
    print(f"Required min size = {int(best_row['required_min_size'])}")
    print(f"Score = {best_row['score']:.4f}")

    return best_k, result_df, valid_df


def compute_soft_labels_in_chunks(X, centers, tau=0.5, chunk_size=50000):
    n = X.shape[0]
    k = centers.shape[0]
    soft_labels = np.zeros((n, k), dtype=np.float32)

    for start in tqdm(range(0, n, chunk_size), desc="Computing soft labels"):
        end = min(start + chunk_size, n)
        X_chunk = X[start:end]

        dist = np.linalg.norm(X_chunk[:, None, :] - centers[None, :, :], axis=-1)
        logits = -dist / tau
        logits = logits - logits.max(axis=1, keepdims=True)

        exp_logits = np.exp(logits)
        soft = exp_logits / exp_logits.sum(axis=1, keepdims=True)
        soft_labels[start:end] = soft.astype(np.float32)

    return soft_labels


def fit_final_kmeans(
    X,
    valid_df,
    best_k,
    out_dir="./proto_k_results",
    random_state=42,
    tau=0.5,
    batch_size=8192,
    soft_chunk_size=50000,
):
    kmeans = MiniBatchKMeans(
        n_clusters=best_k,
        random_state=random_state,
        batch_size=batch_size,
        n_init=20,
        max_iter=500,
        reassignment_ratio=0.01,
        verbose=0,
    )
    labels = kmeans.fit_predict(X)
    centers = kmeans.cluster_centers_.astype(np.float32)
    soft_labels = compute_soft_labels_in_chunks(
        X=X,
        centers=centers,
        tau=tau,
        chunk_size=soft_chunk_size,
    )
    np.save(os.path.join(out_dir, "final_kmeans_centers.npy"), centers)
    np.save(os.path.join(out_dir, "final_kmeans_labels.npy"), labels)
    np.save(os.path.join(out_dir, "final_proto_soft_labels.npy"), soft_labels)
    result_df = valid_df.copy()
    result_df["proto_label"] = labels
    for k in range(best_k):
        result_df[f"proto_prob_{k}"] = soft_labels[:, k]
    result_df.to_csv(os.path.join(out_dir, "final_molecule_proto_labels.csv"), index=False)
    cluster_counts = pd.DataFrame({
        "proto_id": np.arange(best_k),
        "count": np.bincount(labels, minlength=best_k),
    })
    cluster_counts.to_csv(os.path.join(out_dir, "final_cluster_counts.csv"), index=False)
    print("\n========== FINAL CLUSTERING COMPLETED ==========")
    print(f"Final K = {best_k}")
    print("Number of samples in each prototype:")
    print(cluster_counts)
    print(f"\nResults saved to: {out_dir}")
    return labels, centers, soft_labels


def smiles_to_rdkit2d_vector(smi):
    generator = RDKit2DNormalized()

    try:
        vec = np.array(list(generator.process(smi)), dtype=np.float32)
        vec = vec[1:]
        if np.any(np.isnan(vec)) or np.any(np.isinf(vec)):
            raise ValueError(f"The computed SMILES feature vector contains NaN or inf: {smi}")

        return vec

    except Exception as e:
        raise ValueError(f"Failed to compute SMILES features: {smi}, error: {e}")


def compute_two_smiles_similarity(smi1, smi2, metric="cosine"):
    v1 = smiles_to_rdkit2d_vector(smi1)
    v2 = smiles_to_rdkit2d_vector(smi2)

    if v1.shape != v2.shape:
        raise ValueError(f"The feature dimensions of the two molecules do not match: {v1.shape} vs {v2.shape}")

    if metric == "cosine":
        sim = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2) + 1e-12)
        return float(sim)

    elif metric == "euclidean":
        dist = np.linalg.norm(v1 - v2)
        sim = 1.0 / (1.0 + dist)
        return float(sim)

    elif metric == "distance":
        dist = np.linalg.norm(v1 - v2)
        return float(dist)

    else:
        raise ValueError(f"Unsupported metric: {metric}, Available options: cosine / euclidean / distance")

def smiles_list_distance_matrix(smiles_list, metric="cosine"):
    N = len(smiles_list)

    vectors = []
    valid_smiles = []
    valid_indices = []

    for i, smi in enumerate(smiles_list):
        vec = smiles_to_rdkit2d_vector(smi)

        if vec is None:
            print(f"Skipping invalid SMILES: index={i}, smiles={smi}")
            continue

        vectors.append(vec)
        valid_smiles.append(smi)
        valid_indices.append(i)

    if len(vectors) == 0:
        raise ValueError("No SMILES descriptors were computed successfully.")

    vectors = np.vstack(vectors).astype(np.float32)  # [N, D]

    if metric == "cosine":
        norms = np.linalg.norm(vectors, axis=1, keepdims=True)  # [N, 1]
        sim_matrix = vectors @ vectors.T / (norms @ norms.T + 1e-12)
        sim_matrix = np.clip(sim_matrix, -1.0, 1.0)
        dist_matrix = (1.0 - sim_matrix)
        dist_matrix = np.clip(dist_matrix, 0.0, 1.0)

    elif metric == "euclidean":
        dist_matrix = np.zeros((len(vectors), len(vectors)), dtype=np.float32)

        for i in range(len(vectors)):
            diff = vectors - vectors[i]  # [N, D]
            dist = np.linalg.norm(diff, axis=1)
            dist_matrix[i, :] = dist

        # If this matrix will later be used by adaptive_margin_loss, it is recommended to normalize Euclidean distances to [0, 1].
        max_dist = dist_matrix.max()
        if max_dist > 0:
            dist_matrix = dist_matrix / (max_dist + 1e-12)

        dist_matrix = np.clip(dist_matrix, 0.0, 1.0)

    else:
        raise ValueError(f"Unsupported metric: {metric}")

    return dist_matrix.astype(np.float32), valid_smiles, valid_indices


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("--csv_path", type=str, default="./zinc_200w.csv")
    parser.add_argument("--smiles_col", type=str, default="smiles")
    parser.add_argument("--out_dir", type=str, default="./layer1_proto_k_results")
    parser.add_argument("--k_min", type=int, default=40)
    parser.add_argument("--k_max", type=int, default=50)
    parser.add_argument("--random_state", type=int, default=42)
    parser.add_argument("--min_cluster_ratio", type=float, default=0.005)
    parser.add_argument("--min_cluster_size", type=int, default=200)
    parser.add_argument("--batch_size", type=int, default=8192)
    parser.add_argument("--metric_sample_size", type=int, default=50000)
    parser.add_argument("--tau", type=float, default=0.5)
    parser.add_argument("--soft_chunk_size", type=int, default=50000)

    args = parser.parse_args()

    X, valid_df = compute_rdkit2d_vectors(
        csv_path=args.csv_path,
        smiles_col=args.smiles_col,
        out_dir=args.out_dir,
    )

    best_k, result_df, valid_scored_df = explore_kmeans_k(
        X=X,
        k_min=args.k_min,
        k_max=args.k_max,
        out_dir=args.out_dir,
        random_state=args.random_state,
        min_cluster_ratio=args.min_cluster_ratio,
        min_cluster_size=args.min_cluster_size,
        batch_size=args.batch_size,
        metric_sample_size=args.metric_sample_size,
    )

    fit_final_kmeans(
        X=X,
        valid_df=valid_df,
        best_k=best_k,
        out_dir=args.out_dir,
        random_state=args.random_state,
        tau=args.tau,
        batch_size=args.batch_size,
        soft_chunk_size=args.soft_chunk_size,
    )


if __name__ == "__main__":
    main()
