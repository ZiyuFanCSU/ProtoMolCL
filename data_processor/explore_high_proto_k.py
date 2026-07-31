import os
import json
import argparse
import warnings
import numpy as np
import pandas as pd
from tqdm import tqdm
from rdkit import Chem, RDLogger, rdBase
from rdkit.Chem import Descriptors, QED
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import MiniBatchKMeans
from sklearn.metrics import (
    silhouette_score,
    calinski_harabasz_score,
    davies_bouldin_score,
)
RDLogger.DisableLog("rdApp.warning")
rdBase.DisableLog("rdApp.warning")
warnings.filterwarnings("ignore")

SIZE_PROXY_NAMES = [
    "HeavyAtomCount",
]

PHYSICOCHEMICAL_NAMES = [
    "MolLogP",
    "TPSA",
    "QED",
    "MolMR",
]

HBOND_FLEXIBILITY_NAMES = [
    "NumHDonors",
    "NumHAcceptors",
    "NumRotatableBonds",
    "FractionCSP3",
]

RING_SCAFFOLD_NAMES = [
    "RingCount",
    "NumAromaticRings",
    "NumAliphaticRings",
    "NumSaturatedRings",
    "NumAromaticHeterocycles",
    "NumAromaticCarbocycles",
    "NumAliphaticHeterocycles",
    "NumAliphaticCarbocycles",
    "NumSaturatedHeterocycles",
    "NumSaturatedCarbocycles",
]

TOPOLOGY_SHAPE_NAMES = [
    "BalabanJ",
    "HallKierAlpha",
    "Kappa1",
    "Kappa2",
    "Kappa3",
    "Chi0",
    "Chi1",
    "Chi2n",
    "Chi3n",
    "Chi4n",
]

PROPERTY_NAMES = (
    SIZE_PROXY_NAMES
    + PHYSICOCHEMICAL_NAMES
    + HBOND_FLEXIBILITY_NAMES
    + RING_SCAFFOLD_NAMES
    + TOPOLOGY_SHAPE_NAMES
)

FEATURE_GROUPS = {
    "size_proxy": SIZE_PROXY_NAMES,
    "physicochemical": PHYSICOCHEMICAL_NAMES,
    "hbond_flexibility": HBOND_FLEXIBILITY_NAMES,
    "ring_scaffold": RING_SCAFFOLD_NAMES,
    "topology_shape": TOPOLOGY_SHAPE_NAMES,
}


def safe_float(x, default=np.nan):
    try:
        return float(x)
    except Exception:
        return default


def compute_property_vector(smiles):
    """
    Compute balanced global physicochemical-structural descriptors.
    Returns a float32 vector with PROPERTY_NAMES order.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None

    try:
        feats = [
            Descriptors.HeavyAtomCount(mol),
            Descriptors.MolLogP(mol),
            Descriptors.TPSA(mol),
            QED.qed(mol),
            Descriptors.MolMR(mol),
            Descriptors.NumHDonors(mol),
            Descriptors.NumHAcceptors(mol),
            Descriptors.NumRotatableBonds(mol),
            Descriptors.FractionCSP3(mol),

            Descriptors.RingCount(mol),
            Descriptors.NumAromaticRings(mol),
            Descriptors.NumAliphaticRings(mol),
            Descriptors.NumSaturatedRings(mol),
            Descriptors.NumAromaticHeterocycles(mol),
            Descriptors.NumAromaticCarbocycles(mol),
            Descriptors.NumAliphaticHeterocycles(mol),
            Descriptors.NumAliphaticCarbocycles(mol),
            Descriptors.NumSaturatedHeterocycles(mol),
            Descriptors.NumSaturatedCarbocycles(mol),

            Descriptors.BalabanJ(mol),
            Descriptors.HallKierAlpha(mol),
            Descriptors.Kappa1(mol),
            Descriptors.Kappa2(mol),
            Descriptors.Kappa3(mol),
            Descriptors.Chi0(mol),
            Descriptors.Chi1(mol),
            Descriptors.Chi2n(mol),
            Descriptors.Chi3n(mol),
            Descriptors.Chi4n(mol),
        ]
    except Exception:
        return None

    vec = np.array([safe_float(v) for v in feats], dtype=np.float32)
    if np.any(np.isnan(vec)) or np.any(np.isinf(vec)):
        return None
    return vec


def save_feature_metadata(out_dir):
    """Save descriptor names and feature groups for reproducibility."""
    os.makedirs(out_dir, exist_ok=True)
    meta = {
        "feature_names": PROPERTY_NAMES,
        "feature_groups": FEATURE_GROUPS,
        "removed_high_size_redundancy_descriptors": [
            "MolWt",
            "ExactMolWt",
            "LabuteASA",
            "BertzCT",
            "Ipc",
        ],
        "design_note": (
            "Only HeavyAtomCount is kept as a weak size proxy. The remaining features "
            "emphasize physicochemical state, hydrogen bonding, flexibility, ring/scaffold "
            "composition, and topology/shape to reduce mass-dominated clustering."
        ),
    }
    with open(os.path.join(out_dir, "property_feature_metadata.json"), "w", encoding="utf-8") as f:
        json.dump(meta, f, indent=2, ensure_ascii=False)


def compute_property_vectors(csv_path, smiles_col="smiles", out_dir="./property_proto_results"):
    """
    Read SMILES from CSV and compute balanced property-structural descriptor vectors.
    """
    os.makedirs(out_dir, exist_ok=True)
    save_feature_metadata(out_dir)

    df = pd.read_csv(csv_path)
    if smiles_col not in df.columns:
        raise ValueError(f"Column not found in CSV: {smiles_col}. Available columns: {list(df.columns)}")

    smiles_list = df[smiles_col].astype(str).tolist()

    features = []
    valid_indices = []
    failed = []

    for idx, smi in tqdm(list(enumerate(smiles_list)), desc="Computing balanced property-structural descriptors"):
        try:
            vec = compute_property_vector(smi)
            if vec is None:
                failed.append((idx, smi, "invalid_or_descriptor_failed"))
                continue

            features.append(vec)
            valid_indices.append(idx)

        except Exception as e:
            failed.append((idx, smi, str(e)))

    if len(features) == 0:
        raise RuntimeError("No molecular property vectors were computed successfully. Please check the SMILES column and the RDKit environment.")

    X_raw = np.vstack(features).astype(np.float32)

    valid_df = df.iloc[valid_indices].copy()
    valid_df["valid_index"] = np.arange(len(valid_df))

    prop_df = pd.DataFrame(X_raw, columns=PROPERTY_NAMES)
    valid_with_prop_df = pd.concat([valid_df.reset_index(drop=True), prop_df], axis=1)

    np.save(os.path.join(out_dir, "property_features_raw.npy"), X_raw)
    valid_df.to_csv(os.path.join(out_dir, "valid_smiles.csv"), index=False)
    valid_with_prop_df.to_csv(os.path.join(out_dir, "valid_smiles_with_properties.csv"), index=False)

    if len(failed) > 0:
        failed_df = pd.DataFrame(failed, columns=["row_index", "smiles", "reason"])
        failed_df.to_csv(os.path.join(out_dir, "failed_smiles.csv"), index=False)

    print("\n========== Property-structural feature extraction ==========")
    print(f"Number of valid molecules: {len(valid_df)}")
    print(f"Number of failed molecules: {len(failed)}")
    print(f"Feature dimension: {X_raw.shape[1]}")
    print(f"Feature list: {PROPERTY_NAMES}")

    return X_raw, valid_df, valid_with_prop_df


def standardize_features(X_raw, out_dir="./property_proto_results"):
    """
    Standardize continuous descriptors before KMeans.
    """
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X_raw).astype(np.float32)

    np.save(os.path.join(out_dir, "property_features_scaled.npy"), X_scaled)
    np.save(os.path.join(out_dir, "property_scaler_mean.npy"), scaler.mean_)
    np.save(os.path.join(out_dir, "property_scaler_scale.npy"), scaler.scale_)

    print("\n========== Standardization ==========")
    print(f"scaled X shape: {X_scaled.shape}")

    return X_scaled


def save_size_correlation_report(X_raw, out_dir="./property_proto_results", size_feature="HeavyAtomCount"):
    """
    Diagnose whether descriptors are still highly correlated with the retained size proxy.
    This does not remove features automatically; it only saves a report for inspection.
    """
    os.makedirs(out_dir, exist_ok=True)
    df_feat = pd.DataFrame(X_raw, columns=PROPERTY_NAMES)
    if size_feature not in df_feat.columns:
        return None

    corr = df_feat.corr(numeric_only=True)[size_feature].abs().sort_values(ascending=False)
    corr_df = corr.reset_index()
    corr_df.columns = ["feature", f"abs_corr_with_{size_feature}"]
    corr_df["potential_size_shortcut"] = corr_df[f"abs_corr_with_{size_feature}"] >= 0.85
    corr_df.to_csv(os.path.join(out_dir, "feature_size_correlation_report.csv"), index=False)

    print("\n========== Size-shortcut diagnostic ==========")
    print(corr_df.head(12))
    return corr_df


def sample_for_metric(X, labels, sample_size=20000, random_state=42):
    """
    Sample data for clustering metrics to avoid expensive full silhouette on large datasets.
    """
    n = X.shape[0]
    if n <= sample_size:
        return X, labels

    rng = np.random.default_rng(random_state)
    idx = rng.choice(n, size=sample_size, replace=False)
    return X[idx], labels[idx]


def explore_kmeans_k(
    X,
    k_min=5,
    k_max=30,
    out_dir="./property_proto_results",
    random_state=42,
    min_cluster_ratio=0.005,
    min_cluster_size=200,
    batch_size=8192,
    metric_sample_size=20000,
):
    """
    Search K for Layer-5 balanced property-structural prototypes with MiniBatchKMeans.
    Metrics are computed on sampled data for scalability.
    """
    os.makedirs(out_dir, exist_ok=True)

    results = []
    n_samples = X.shape[0]

    for k in tqdm(range(k_min, k_max + 1), desc="Exploring K"):
        if k >= n_samples:
            continue
        if k < 2:
            warnings.warn(f"K={k} was skipped because the Silhouette, CH, and DB metrics require at least two clusters.")
            continue

        try:
            kmeans = MiniBatchKMeans(
                n_clusters=k,
                random_state=random_state,
                batch_size=batch_size,
                n_init=10,
                max_iter=300,
                reassignment_ratio=0.01,
            )

            labels = kmeans.fit_predict(X)
            counts = np.bincount(labels, minlength=k)

            required_min_size = max(min_cluster_size, int(n_samples * min_cluster_ratio))
            valid_cluster_size = counts.min() >= required_min_size

            X_metric, labels_metric = sample_for_metric(
                X,
                labels,
                sample_size=metric_sample_size,
                random_state=random_state,
            )

            sil = silhouette_score(X_metric, labels_metric)
            ch = calinski_harabasz_score(X_metric, labels_metric)
            db = davies_bouldin_score(X_metric, labels_metric)

            results.append({
                "K": k,
                "silhouette_sample": sil,
                "calinski_harabasz_sample": ch,
                "davies_bouldin_sample": db,
                "inertia": kmeans.inertia_,
                "min_cluster_size": int(counts.min()),
                "max_cluster_size": int(counts.max()),
                "required_min_size": int(required_min_size),
                "valid_cluster_size": bool(valid_cluster_size),
                "cluster_counts": counts.tolist(),
            })

        except Exception as e:
            warnings.warn(f"Clustering failed for K={k}: {e}")

    result_df = pd.DataFrame(results)
    if len(result_df) == 0:
        raise RuntimeError("No clustering run completed successfully for any K. Please check the input features or increase k_max.")

    result_df.to_csv(os.path.join(out_dir, "property_k_search_all.csv"), index=False)

    valid_df = result_df[result_df["valid_cluster_size"] == True].copy()
    if len(valid_df) == 0:
        print("\nWarning: No K satisfies the minimum cluster-size requirement. Selection will be made from all K values.")
        valid_df = result_df.copy()

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

    valid_df.to_csv(os.path.join(out_dir, "property_k_search_valid_scored.csv"), index=False)

    best_row = valid_df.sort_values("score", ascending=False).iloc[0]
    best_k = int(best_row["K"])

    print("\n========== K Search Results ==========")
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

    print("\n========== Recommended K ==========")
    print(f"Best K = {best_k}")
    print(f"Silhouette(sample) = {best_row['silhouette_sample']:.4f}")
    print(f"Calinski-Harabasz(sample) = {best_row['calinski_harabasz_sample']:.4f}")
    print(f"Davies-Bouldin(sample) = {best_row['davies_bouldin_sample']:.4f}")
    print(f"Min cluster size = {int(best_row['min_cluster_size'])}")
    print(f"Required min size = {int(best_row['required_min_size'])}")
    print(f"Score = {best_row['score']:.4f}")

    return best_k, result_df, valid_df


def compute_soft_labels_chunked(X, centers, tau=0.5, chunk_size=50000):
    """
    Compute soft labels by chunks to avoid memory peaks.
    Uses squared Euclidean distance with matrix multiplication for speed.
    """
    if tau <= 0:
        raise ValueError("--tau must be positive.")

    n = X.shape[0]
    k = centers.shape[0]
    soft_labels = np.zeros((n, k), dtype=np.float32)

    centers = centers.astype(np.float32)
    center_norm = np.sum(centers ** 2, axis=1, keepdims=True).T  # [1, K]

    for start in tqdm(range(0, n, chunk_size), desc="Computing soft labels"):
        end = min(start + chunk_size, n)
        X_chunk = X[start:end].astype(np.float32)
        x_norm = np.sum(X_chunk ** 2, axis=1, keepdims=True)  # [B, 1]
        dist2 = x_norm + center_norm - 2.0 * np.dot(X_chunk, centers.T)
        dist2 = np.maximum(dist2, 0.0)
        logits = -dist2 / tau
        logits = logits - logits.max(axis=1, keepdims=True)
        exp_logits = np.exp(logits)
        soft = exp_logits / (exp_logits.sum(axis=1, keepdims=True) + 1e-12)
        soft_labels[start:end] = soft.astype(np.float32)

    return soft_labels


def fit_final_kmeans(
    X,
    valid_df,
    best_k,
    out_dir="./property_proto_results",
    random_state=42,
    tau=0.5,
    batch_size=8192,
    soft_chunk_size=50000,
    save_soft_csv=False,
    soft_csv_sample_size=20000,
):
    """
    Fit final MiniBatchKMeans and save hard labels, centers and soft labels.
    By default, only .npy soft labels are saved to avoid huge CSV files.
    """
    kmeans = MiniBatchKMeans(
        n_clusters=best_k,
        random_state=random_state,
        batch_size=batch_size,
        n_init=20,
        max_iter=500,
        reassignment_ratio=0.01,
    )

    labels = kmeans.fit_predict(X)
    centers = kmeans.cluster_centers_.astype(np.float32)

    soft_labels = compute_soft_labels_chunked(
        X,
        centers,
        tau=tau,
        chunk_size=soft_chunk_size,
    )

    np.save(os.path.join(out_dir, "final_property_structural_proto_centers.npy"), centers)
    np.save(os.path.join(out_dir, "final_property_structural_proto_labels.npy"), labels)
    np.save(os.path.join(out_dir, "final_property_structural_proto_soft_labels.npy"), soft_labels)

    # Save compact hard-label CSV for inspection and downstream alignment.
    hard_df = valid_df.copy()
    hard_df["property_structural_proto_label"] = labels
    hard_df.to_csv(os.path.join(out_dir, "final_molecule_property_structural_proto_labels.csv"), index=False)

    # Optional sampled soft-label CSV for inspection only. Full soft labels are already saved as .npy.
    if save_soft_csv:
        sample_n = min(int(soft_csv_sample_size), len(valid_df))
        rng = np.random.default_rng(random_state)
        sample_idx = np.sort(rng.choice(len(valid_df), size=sample_n, replace=False))
        soft_df = valid_df.iloc[sample_idx].copy()
        soft_df["property_structural_proto_label"] = labels[sample_idx]
        for k in range(best_k):
            soft_df[f"property_structural_proto_prob_{k}"] = soft_labels[sample_idx, k]
        soft_df.to_csv(os.path.join(out_dir, "sampled_molecule_property_structural_proto_soft_labels.csv"), index=False)

    cluster_counts = pd.DataFrame({
        "property_structural_proto_id": np.arange(best_k),
        "count": np.bincount(labels, minlength=best_k),
    })
    cluster_counts.to_csv(os.path.join(out_dir, "final_property_structural_proto_counts.csv"), index=False)

    with open(os.path.join(out_dir, "final_output_info.json"), "w", encoding="utf-8") as f:
        json.dump({
            "hard_label_csv": "final_molecule_property_structural_proto_labels.csv",
            "centers_npy": "final_property_structural_proto_centers.npy",
            "hard_labels_npy": "final_property_structural_proto_labels.npy",
            "soft_labels_npy": "final_property_structural_proto_soft_labels.npy",
            "soft_label_shape": [int(soft_labels.shape[0]), int(soft_labels.shape[1])],
            "tau": float(tau),
            "save_soft_csv": bool(save_soft_csv),
            "soft_csv_note": "CSV is sampled only when --save_soft_csv is used; full soft labels are stored in .npy.",
        }, f, indent=2, ensure_ascii=False)

    print("\n========== Final Layer-5 Property-Structural Prototype Clustering Completed ==========")
    print(f"Final K = {best_k}")
    print("Number of samples in each prototype:")
    print(cluster_counts)
    print(f"\nResults have been saved to: {out_dir}")

    return labels, centers, soft_labels


from sklearn.preprocessing import StandardScaler
from sklearn.metrics.pairwise import cosine_similarity
import numpy as np


def compute_smiles_distance_matrix(smiles_list, standardize=True):
    """
    Inputs:
        smiles_list: list[str], a list of SMILES strings
        standardize: whether to standardize the descriptors; defaults to True, consistent with the preprocessing before KMeans in this script

    Outputs:
        dist_matrix: numpy.ndarray, shape = [N, N]
                     Distance matrix, where dist_matrix[i, j] is the distance between molecules i and j
                     0 indicates the highest similarity, and 1 indicates the lowest similarity
        valid_smiles: SMILES strings whose descriptors were computed successfully
        valid_indices: indices of successfully processed SMILES strings in the original list
    """

    features = []
    valid_smiles = []
    valid_indices = []

    for i, smi in enumerate(smiles_list):
        vec = compute_property_vector(smi)

        if vec is None:
            print(f"Skipping invalid SMILES: index={i}, smiles={smi}")
            continue

        features.append(vec)
        valid_smiles.append(smi)
        valid_indices.append(i)

    if len(features) == 0:
        raise ValueError("No SMILES string was successfully converted into descriptors.")

    X = np.vstack(features).astype(np.float32)

    # Consistent with the original code: standardize first, then compute distances.
    if standardize:
        scaler = StandardScaler()
        X = scaler.fit_transform(X).astype(np.float32)

    # Cosine similarity theoretically ranges from -1 to 1.
    sim_matrix = cosine_similarity(X)

    # Guard against floating-point errors.
    sim_matrix = np.clip(sim_matrix, -1.0, 1.0)

    # Convert similarity to a distance in the range [0, 1]:
    # sim = 1  -> dist = 0, highest similarity
    # sim = -1 -> dist = 1, lowest similarity
    dist_matrix = (1.0 - sim_matrix) / 2.0

    # Clip again to ensure the result remains within [0, 1].
    dist_matrix = np.clip(dist_matrix, 0.0, 1.0)

    # Force the diagonal entries to zero.
    np.fill_diagonal(dist_matrix, 0.0)

    return dist_matrix.astype(np.float32), valid_smiles, valid_indices

def main():
    parser = argparse.ArgumentParser(
        description="Balanced global physicochemical-structural prototype construction."
    )

    parser.add_argument("--csv_path", type=str, default="./zinc_200w.csv")
    parser.add_argument("--smiles_col", type=str, default="smiles")
    parser.add_argument("--out_dir", type=str, default="./layer5_proto_k_results")
    parser.add_argument("--k_min", type=int, default=1)
    parser.add_argument("--k_max", type=int, default=6)
    parser.add_argument("--random_state", type=int, default=42)
    parser.add_argument("--min_cluster_ratio", type=float, default=0.005)
    parser.add_argument("--min_cluster_size", type=int, default=200)
    parser.add_argument("--batch_size", type=int, default=8192)
    parser.add_argument("--metric_sample_size", type=int, default=20000)
    parser.add_argument("--soft_chunk_size", type=int, default=50000)
    parser.add_argument("--tau", type=float, default=0.5)
    parser.add_argument("--save_soft_csv", action="store_true", help="Save a sampled soft-label CSV for inspection. Full soft labels are always saved as .npy.")
    parser.add_argument("--soft_csv_sample_size", type=int, default=20000)

    args = parser.parse_args()

    X_raw, valid_df, _ = compute_property_vectors(
        csv_path=args.csv_path,
        smiles_col=args.smiles_col,
        out_dir=args.out_dir,
    )

    save_size_correlation_report(
        X_raw=X_raw,
        out_dir=args.out_dir,
        size_feature="HeavyAtomCount",
    )

    X_scaled = standardize_features(
        X_raw=X_raw,
        out_dir=args.out_dir,
    )

    best_k, _, _ = explore_kmeans_k(
        X=X_scaled,
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
        X=X_scaled,
        valid_df=valid_df,
        best_k=best_k,
        out_dir=args.out_dir,
        random_state=args.random_state,
        tau=args.tau,
        batch_size=args.batch_size,
        soft_chunk_size=args.soft_chunk_size,
        save_soft_csv=args.save_soft_csv,
        soft_csv_sample_size=args.soft_csv_sample_size,
    )


if __name__ == "__main__":
    main()
