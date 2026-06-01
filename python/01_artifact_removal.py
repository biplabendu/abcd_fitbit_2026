# %% [markdown]
# # 01 — Artifact Removal
#
# Removes physiologically implausible values from each configured signal.
# Signals, their artifact bounds, and spike thresholds are defined in
# `config.SIGNALS`. Optionally restricts to a curated participant ID list.
#
# **Input:**  `dev/data/fitbit-summaries/activity_120m.parquet`
# **Output:** `data/processed/clean/clean_data.parquet`

# %% Imports and setup
import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from config import (
    RAW_PARQUET, CLEAN_PARQUET, LOGS_DIR,
    SESSIONS,
    APPLY_ID_FILTER, ID_FILTER_FILE, ID_FILTER_COLUMN,
    SIGNALS, SIGNAL_NAMES, VALIDITY_SIGNAL,
    COL_ID, COL_SESSION, COL_DAY, COL_START, COL_WKND,
)

LOGS_DIR.mkdir(parents=True, exist_ok=True)
CLEAN_PARQUET.parent.mkdir(parents=True, exist_ok=True)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)s  %(message)s",
    handlers=[
        logging.FileHandler(LOGS_DIR / "01_artifact_removal.log"),
        logging.StreamHandler(),
    ],
)
log = logging.getLogger(__name__)

# %% Load raw data
log.info(f"Loading {RAW_PARQUET}")
raw = pd.read_parquet(RAW_PARQUET)
log.info(f"  Raw rows: {len(raw):,}  ({raw[COL_ID].nunique():,} participants)")

df = raw[raw[COL_SESSION].isin(SESSIONS)].copy()
log.info(f"  After session filter {SESSIONS}: {len(df):,} rows, "
         f"{df[COL_ID].nunique():,} participants")

# %% Restrict to curated participant IDs
if APPLY_ID_FILTER:
    keep_ids = pd.read_csv(ID_FILTER_FILE)[ID_FILTER_COLUMN].astype(str).unique()
    df = df[df[COL_ID].isin(keep_ids)].copy()
    log.info(f"  After ID filter ({ID_FILTER_FILE.name}): {len(df):,} rows, "
             f"{df[COL_ID].nunique():,} of {len(keep_ids):,} listed participants")
else:
    log.info("  ID filter disabled — running on all participants")

# %% Parse timestamps and derive start_hour
df[COL_START] = pd.to_datetime(df[COL_START])
df["start_hour"] = df[COL_START].dt.hour   # 0, 2, 4, …, 22

# %% Artifact removal — per signal
for sig, meta in SIGNALS.items():
    n_before = df[sig].notna().sum()

    # Hard range bounds
    lo, hi = meta.get("clip_min"), meta.get("clip_max")
    if lo is not None:
        df.loc[df[sig] < lo, sig] = np.nan
    if hi is not None:
        df.loc[df[sig] > hi, sig] = np.nan

    # Within participant-day spike removal (activity-type signals)
    spike_sd = meta.get("spike_sd")
    if spike_sd is not None:
        def _spike_mask(grp: pd.DataFrame, _sig=sig, _sd=spike_sd) -> pd.Series:
            mu, sd = grp[_sig].mean(), grp[_sig].std()
            if pd.isna(sd) or sd == 0:
                return pd.Series(False, index=grp.index)
            return grp[_sig] > mu + _sd * sd

        spikes = (
            df.groupby([COL_ID, COL_SESSION, COL_DAY], group_keys=False)
            .apply(_spike_mask)
        )
        df.loc[spikes, sig] = np.nan

    n_removed = n_before - df[sig].notna().sum()
    log.info(f"  {sig}: removed {n_removed:,} "
             f"({n_removed / max(n_before, 1):.1%})  [{meta['kind']}]")

# TODO: sleep_stages — add Fitbit sleep stage code validation here once raw
#       per-minute stage data is available (currently only min_slp is used).

# %% Save
cols = ([COL_ID, COL_SESSION, COL_DAY, COL_START, "start_hour", COL_WKND]
        + SIGNAL_NAMES
        + ([VALIDITY_SIGNAL] if VALIDITY_SIGNAL not in SIGNAL_NAMES else []))
out = df[cols]
out.to_parquet(CLEAN_PARQUET, index=False)
log.info(f"Saved → {CLEAN_PARQUET}  ({len(out):,} rows, signals={SIGNAL_NAMES})")
