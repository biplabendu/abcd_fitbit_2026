# %% [markdown]
# # 05 — Dimensionality Reduction
#
# PCA (retain 95 % variance) → UMAP per session.
# Produces a 10-component embedding for clustering and a
# 2-component embedding for visualisation.
#
# **Input:**  `data/handoff/features_{yr2,yr6}.csv`
# **Output:** `data/processed/umap_models/` (PCA + UMAP models, embeddings)
#             `data/handoff/umap_coords.csv`

# %% Imports and setup
import logging
import pickle
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
import umap

sys.path.insert(0, str(Path(__file__).parent))
from config import (
    HANDOFF_DIR, MODELS_DIR, LOGS_DIR,
    SESSIONS, SESSION_LABELS,
    COL_ID, COL_SESSION,
    PCA_VARIANCE_THRESHOLD,
    UMAP_N_NEIGHBORS_FINAL, UMAP_MIN_DIST_FINAL,
    RANDOM_STATE,
)

LOGS_DIR.mkdir(parents=True, exist_ok=True)
MODELS_DIR.mkdir(parents=True, exist_ok=True)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)s  %(message)s",
    handlers=[
        logging.FileHandler(LOGS_DIR / "05_dimensionality_reduction.log"),
        logging.StreamHandler(),
    ],
)
log = logging.getLogger(__name__)

# %% Load features
all_dfs = []
for sess in SESSIONS:
    label = SESSION_LABELS[sess]
    path  = HANDOFF_DIR / f"features_{label}.csv"
    wave  = pd.read_csv(path)
    assert COL_ID in wave.columns, f"Missing {COL_ID} in {path}"
    all_dfs.append(wave)
    log.info(f"  Loaded {label}: {len(wave):,} rows")

features = pd.concat(all_dfs, ignore_index=True)
z_cols   = [c for c in features.columns if c.startswith("z_")]
log.info(f"  Z-scored features: {len(z_cols)}")

# %% PCA + UMAP per session
umap_rows = []

for sess in SESSIONS:
    label = SESSION_LABELS[sess]
    mask  = features[COL_SESSION] == sess
    sub   = features[mask].copy()

    X = sub[z_cols].values

    # Drop rows with any NaN remaining in the feature matrix
    valid_rows = ~np.isnan(X).any(axis=1)
    n_dropped  = int((~valid_rows).sum())
    if n_dropped:
        log.warning(f"  {label}: dropping {n_dropped} rows with NaN features")
    X     = X[valid_rows]
    ids   = sub[valid_rows][COL_ID].values

    # -- PCA ----------------------------------------------------------------
    pca    = PCA(n_components=PCA_VARIANCE_THRESHOLD, random_state=RANDOM_STATE)
    X_pca  = pca.fit_transform(X)
    log.info(f"  {label}: PCA {X.shape[1]}→{X_pca.shape[1]} components "
             f"({pca.explained_variance_ratio_.sum():.1%} variance)")

    # -- UMAP 10-D (for clustering) -----------------------------------------
    u10 = umap.UMAP(
        n_components=10,
        n_neighbors=UMAP_N_NEIGHBORS_FINAL,
        min_dist=UMAP_MIN_DIST_FINAL,
        random_state=RANDOM_STATE,
        metric="euclidean",
    )
    X_umap10 = u10.fit_transform(X_pca)

    # -- UMAP 2-D (for visualisation only) ----------------------------------
    u2 = umap.UMAP(
        n_components=2,
        n_neighbors=UMAP_N_NEIGHBORS_FINAL,
        min_dist=UMAP_MIN_DIST_FINAL,
        random_state=RANDOM_STATE,
        metric="euclidean",
    )
    X_umap2 = u2.fit_transform(X_pca)
    log.info(f"  {label}: UMAP done")

    # -- Persist models and embeddings --------------------------------------
    with open(MODELS_DIR / f"pca_{label}.pkl",    "wb") as f: pickle.dump(pca, f)
    with open(MODELS_DIR / f"umap10_{label}.pkl", "wb") as f: pickle.dump(u10, f)
    with open(MODELS_DIR / f"umap2_{label}.pkl",  "wb") as f: pickle.dump(u2,  f)

    np.save(MODELS_DIR / f"umap10_embedding_{label}.npy", X_umap10)
    np.save(MODELS_DIR / f"ids_{label}.npy",              ids)

    # -- Collect 2-D coords for handoff ------------------------------------
    for j, pid in enumerate(ids):
        umap_rows.append({
            COL_ID: pid, COL_SESSION: sess, "wave_label": label,
            "umap_1": float(X_umap2[j, 0]),
            "umap_2": float(X_umap2[j, 1]),
        })

# %% Save UMAP visualisation coordinates
umap_df = pd.DataFrame(umap_rows)
umap_df.to_csv(HANDOFF_DIR / "umap_coords.csv", index=False)
log.info(f"Saved UMAP coords → {HANDOFF_DIR / 'umap_coords.csv'}")
