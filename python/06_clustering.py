# %% [markdown]
# # 06 — Clustering, Stability, and Cross-Wave Alignment
#
# **Primary:**   HDBSCAN (parameter sweep, best by silhouette; noise resolved by 1-NN)
# **Secondary:** k-medoids (k = 3–7, diagnostic only — not used for exported labels)
# **Soft:**      GMM with k = number of HDBSCAN clusters
# **Stability:** Bootstrap ARI (100 iterations, 80 % resample)
# **Alignment:** Hungarian algorithm on feature-space centroid cosine distance
#
# **Input:**  `data/processed/umap_models/umap10_embedding_{label}.npy`
#             `data/handoff/features_{yr2,yr6}.csv`
# **Output:** `data/handoff/cluster_assignments.csv`
#             `data/handoff/gmm_probabilities.csv`
#             `data/handoff/cluster_centroids.csv`

# %% Imports and setup
import logging
import pickle
import sys
from itertools import product
from pathlib import Path

import numpy as np
import pandas as pd
from hdbscan import HDBSCAN
from scipy.optimize import linear_sum_assignment
from scipy.spatial.distance import cdist
from sklearn.metrics import (
    adjusted_rand_score,
    calinski_harabasz_score,
    davies_bouldin_score,
    silhouette_score,
)
from sklearn.mixture import GaussianMixture
from sklearn.neighbors import NearestNeighbors
try:
    from sklearn_extra.cluster import KMedoids
    KMEDOIDS_PRECOMPUTED = False
except ImportError:
    from kmedoids import KMedoids
    KMEDOIDS_PRECOMPUTED = True  # kmedoids package needs a distance matrix, not raw features

sys.path.insert(0, str(Path(__file__).parent))
from config import (
    HANDOFF_DIR, MODELS_DIR, LOGS_DIR,
    SESSIONS, SESSION_LABELS,
    COL_ID, COL_SESSION,
    HDBSCAN_MIN_CLUSTER_SIZES, HDBSCAN_MIN_SAMPLES,
    KMEDOIDS_K_RANGE,
    N_BOOTSTRAP, BOOTSTRAP_FRACTION,
    RANDOM_STATE,
)

LOGS_DIR.mkdir(parents=True, exist_ok=True)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)s  %(message)s",
    handlers=[
        logging.FileHandler(LOGS_DIR / "06_clustering.log"),
        logging.StreamHandler(),
    ],
)
log = logging.getLogger(__name__)


def _kmedoids_fit(X: np.ndarray, k: int, random_state: int) -> np.ndarray:
    """Fit k-medoids, handling both sklearn_extra (raw features) and
    kmedoids package (requires precomputed distance matrix) APIs."""
    if KMEDOIDS_PRECOMPUTED:
        D = cdist(X, X, metric="euclidean").astype(np.float32)
        return KMedoids(n_clusters=k, random_state=random_state).fit_predict(D)
    return KMedoids(n_clusters=k, random_state=random_state).fit_predict(X)


# %% Load embeddings and feature matrices
embeddings     = {}
ids_by_sess    = {}
features_by_sess = {}

for sess in SESSIONS:
    label = SESSION_LABELS[sess]
    embeddings[sess]      = np.load(MODELS_DIR / f"umap10_embedding_{label}.npy")
    ids_by_sess[sess]     = np.load(MODELS_DIR / f"ids_{label}.npy", allow_pickle=True)
    features_by_sess[sess] = pd.read_csv(HANDOFF_DIR / f"features_{label}.csv")
    log.info(f"  {label}: {embeddings[sess].shape[0]} participants, "
             f"{embeddings[sess].shape[1]} UMAP dims")

# %% HDBSCAN parameter sweep
log.info("HDBSCAN sweep...")
hdbscan_best = {}

for sess in SESSIONS:
    label = SESSION_LABELS[sess]
    X     = embeddings[sess]
    best  = {"sil": -np.inf, "labels": None, "params": None}

    for mcs, ms in product(HDBSCAN_MIN_CLUSTER_SIZES, HDBSCAN_MIN_SAMPLES):
        labels = HDBSCAN(min_cluster_size=mcs, min_samples=ms).fit_predict(X)
        mask   = labels != -1
        k      = len(set(labels[mask]))
        if k < 2 or mask.sum() < 10:
            continue
        sil = silhouette_score(X[mask], labels[mask])
        log.info(f"  {label} mcs={mcs} ms={ms}: k={k} "
                 f"noise={1 - mask.mean():.1%} sil={sil:.3f}")
        if sil > best["sil"]:
            best = {"sil": sil, "labels": labels, "params": (mcs, ms)}

    hdbscan_best[sess] = best
    log.info(f"  {label} best: {best['params']}  sil={best['sil']:.3f}")

# %% Resolve HDBSCAN noise labels → primary cluster labels (1-NN assignment)
log.info("Resolving HDBSCAN noise labels...")
primary_labels = {}
for sess in SESSIONS:
    label = SESSION_LABELS[sess]
    raw   = hdbscan_best[sess]["labels"].copy()
    noise = raw == -1
    if noise.any() and (~noise).sum() >= 2:
        X  = embeddings[sess]
        nn = NearestNeighbors(n_neighbors=1).fit(X[~noise])
        _, idx = nn.kneighbors(X[noise])
        raw[noise] = raw[~noise][idx.flatten()]
        log.info(f"  {label}: resolved {int(noise.sum())} noise points")
    else:
        log.info(f"  {label}: 0 noise points")
    primary_labels[sess] = raw

# %% k-medoids sweep (diagnostic — metrics only, labels not used for exports)
log.info("k-medoids sweep...")
kmedoids_results = {}

for sess in SESSIONS:
    label   = SESSION_LABELS[sess]
    X       = embeddings[sess]
    metrics = []

    for k in KMEDOIDS_K_RANGE:
        labels = _kmedoids_fit(X, k, RANDOM_STATE)
        metrics.append({
            "k":                  k,
            "silhouette":         silhouette_score(X, labels),
            "davies_bouldin":     davies_bouldin_score(X, labels),
            "calinski_harabasz":  calinski_harabasz_score(X, labels),
            "labels":             labels,
        })
        log.info(f"  {label} k={k}: sil={metrics[-1]['silhouette']:.3f} "
                 f"db={metrics[-1]['davies_bouldin']:.3f}")

    best_idx             = int(np.argmax([m["silhouette"] for m in metrics]))
    kmedoids_results[sess] = {"metrics": metrics, "best": metrics[best_idx]}
    log.info(f"  {label} best k={kmedoids_results[sess]['best']['k']}")

# %% GMM (soft assignments using k from primary HDBSCAN labels)
log.info("Fitting GMM...")
gmm_rows = {}

for sess in SESSIONS:
    label = SESSION_LABELS[sess]
    X     = embeddings[sess]
    ids   = ids_by_sess[sess]
    k     = len(np.unique(primary_labels[sess]))

    gmm   = GaussianMixture(n_components=k, random_state=RANDOM_STATE, n_init=5)
    gmm.fit(X)
    probs = gmm.predict_proba(X)   # (n, k)
    hard  = gmm.predict(X)

    gmm_rows[sess] = []
    for j, pid in enumerate(ids):
        row = {COL_ID: pid, COL_SESSION: sess, "wave_label": label}
        for ki in range(k):
            row[f"prob_cluster_{ki + 1}"] = float(probs[j, ki])
        gmm_rows[sess].append(row)

    log.info(f"  {label} GMM k={k} done")

# %% Bootstrap stability (ARI)
log.info(f"Bootstrap stability ({N_BOOTSTRAP} iterations)...")

for sess in SESSIONS:
    label = SESSION_LABELS[sess]
    X     = embeddings[sess]
    k     = len(np.unique(primary_labels[sess]))
    full  = primary_labels[sess]
    n     = X.shape[0]
    rng   = np.random.default_rng(RANDOM_STATE)
    aris  = []

    for b in range(N_BOOTSTRAP):
        idx  = rng.choice(n, size=int(n * BOOTSTRAP_FRACTION), replace=False)
        boot = _kmedoids_fit(X[idx], k, int(b))
        aris.append(adjusted_rand_score(full[idx], boot))

    log.info(f"  {label} ARI: {np.mean(aris):.3f} ± {np.std(aris):.3f}")

# %% Cross-wave cluster alignment (Hungarian algorithm)
log.info("Cross-wave alignment...")

z_cols_ref     = None
centroid_rows  = []
labels_raw     = {}
labels_aligned = {}

for sess in SESSIONS:
    label  = SESSION_LABELS[sess]
    feat   = features_by_sess[sess]
    ids    = ids_by_sess[sess]
    lbl    = primary_labels[sess]
    z_cols = [c for c in feat.columns if c.startswith("z_")]

    if z_cols_ref is None:
        z_cols_ref = z_cols

    labels_raw[sess] = lbl
    X_feat = feat.set_index(COL_ID).reindex(ids)[z_cols].fillna(0).values

    for k_val in np.unique(lbl):
        centroid = X_feat[lbl == k_val].mean(axis=0)
        row = {"cluster": int(k_val), COL_SESSION: sess, "wave_label": label}
        for j, col in enumerate(z_cols):
            row[col] = float(centroid[j])
        centroid_rows.append(row)

# Align all sessions to the first session's labelling
ref_sess = SESSIONS[0]

def _centroids(feat, ids, labels, z_cols):
    X = feat.set_index(COL_ID).reindex(ids)[z_cols].fillna(0).values
    return np.array([X[labels == k].mean(axis=0) for k in sorted(np.unique(labels))])

ref_c = _centroids(features_by_sess[ref_sess], ids_by_sess[ref_sess],
                    labels_raw[ref_sess], z_cols_ref)
labels_aligned[ref_sess] = labels_raw[ref_sess].copy()

alignment_rows = []   # cross-wave centroid similarity for each matched pair
for sess in SESSIONS[1:]:
    tgt_c = _centroids(features_by_sess[sess], ids_by_sess[sess],
                        labels_raw[sess], z_cols_ref)
    cos   = cdist(ref_c, tgt_c, metric="cosine")
    ri, ti = linear_sum_assignment(cos)
    remap  = {int(old): int(new) for new, old in zip(ri, ti)}
    labels_aligned[sess] = np.array([remap.get(int(l), int(l))
                                      for l in labels_raw[sess]])
    log.info(f"  Alignment {SESSION_LABELS[ref_sess]}→{SESSION_LABELS[sess]}: {remap}")

    # Record matched-pair cosine similarity (1 - cosine distance). High values
    # justify pooling the aligned label across waves; low values mean the
    # phenotype is not preserved and cross-wave cluster effects are not
    # interpretable (see ANALYSIS_PLAN §3.4 / Supp Table 3).
    for ref_k, tgt_k in zip(ri, ti):
        alignment_rows.append({
            "ref_session":     ref_sess,
            "ref_wave":        SESSION_LABELS[ref_sess],
            "target_session":  sess,
            "target_wave":     SESSION_LABELS[sess],
            "ref_cluster":     int(ref_k),
            "target_cluster_raw":     int(tgt_k),
            "aligned_cluster": int(ref_k),
            "cosine_similarity": float(1.0 - cos[ref_k, tgt_k]),
        })

# %% Save outputs
assign_rows = []
for sess in SESSIONS:
    label = SESSION_LABELS[sess]
    ids   = ids_by_sess[sess]
    for j, pid in enumerate(ids):
        assign_rows.append({
            COL_ID:            pid,
            COL_SESSION:       sess,
            "wave_label":      label,
            "cluster_hard":    int(labels_raw[sess][j]),
            "cluster_aligned": int(labels_aligned[sess][j]),
        })

pd.DataFrame(assign_rows).to_csv(
    HANDOFF_DIR / "cluster_assignments.csv", index=False)
log.info(f"Saved cluster assignments → {HANDOFF_DIR / 'cluster_assignments.csv'}")

gmm_all = [row for rows in gmm_rows.values() for row in rows]
pd.DataFrame(gmm_all).to_csv(
    HANDOFF_DIR / "gmm_probabilities.csv", index=False)
log.info(f"Saved GMM probabilities → {HANDOFF_DIR / 'gmm_probabilities.csv'}")

pd.DataFrame(centroid_rows).to_csv(
    HANDOFF_DIR / "cluster_centroids.csv", index=False)
log.info(f"Saved cluster centroids → {HANDOFF_DIR / 'cluster_centroids.csv'}")

pd.DataFrame(alignment_rows).to_csv(
    HANDOFF_DIR / "cluster_alignment.csv", index=False)
log.info(f"Saved cross-wave alignment → {HANDOFF_DIR / 'cluster_alignment.csv'}")

log.info("Clustering complete.")
