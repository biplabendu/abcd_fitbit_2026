# %% [markdown]
# # 01 — Artifact Removal
#
# Removes physiologically implausible values from HR, steps, and sleep signals.
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
    COL_ID, COL_SESSION, COL_DAY, COL_START, COL_WKND,
    COL_HR, COL_SLEEP, COL_STEPS,
    HR_MIN, HR_MAX, STEPS_SPIKE_SD, SLEEP_MIN, SLEEP_MAX,
)
from utils.wear_quality import flag_nonwear

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
log.info(f"  Raw rows: {len(raw):,}")

df = raw[raw[COL_SESSION].isin(SESSIONS)].copy()
log.info(f"  After session filter {SESSIONS}: {len(df):,} rows, "
         f"{df[COL_ID].nunique():,} participants")

# %% Parse timestamps and derive start_hour
df[COL_START] = pd.to_datetime(df[COL_START])
df["start_hour"] = df[COL_START].dt.hour   # 0, 2, 4, …, 22

# %% Artifact removal — HR
n_before = df[COL_HR].notna().sum()

out_of_range = (df[COL_HR] < HR_MIN) | (df[COL_HR] > HR_MAX)
df.loc[out_of_range, COL_HR] = np.nan

nonwear = flag_nonwear(df)
df.loc[nonwear, COL_HR] = np.nan
# TODO: sleep_stages — non-wear detection is a stub (0 slots flagged).
#       Implement once raw per-minute stage data is available.

n_removed = n_before - df[COL_HR].notna().sum()
log.info(f"  HR removed: {n_removed:,} ({n_removed / max(n_before, 1):.1%})")

# %% Artifact removal — Steps
n_before = df[COL_STEPS].notna().sum()

neg = df[COL_STEPS] < 0


def _spike_mask(grp: pd.DataFrame) -> pd.Series:
    mu = grp[COL_STEPS].mean()
    sd = grp[COL_STEPS].std()
    if pd.isna(sd) or sd == 0:
        return pd.Series(False, index=grp.index)
    return grp[COL_STEPS] > mu + STEPS_SPIKE_SD * sd


spikes = (
    df.groupby([COL_ID, COL_SESSION, COL_DAY], group_keys=False)
    .apply(_spike_mask)
)
df.loc[neg | spikes, COL_STEPS] = np.nan

n_removed = n_before - df[COL_STEPS].notna().sum()
log.info(f"  Steps removed: {n_removed:,} ({n_removed / max(n_before, 1):.1%})")

# %% Artifact removal — Sleep (min_slp range check)
# TODO: sleep_stages — add Fitbit sleep stage code validation here once
#       raw per-minute data is available.
n_before = df[COL_SLEEP].notna().sum()

out_of_range = (df[COL_SLEEP] < SLEEP_MIN) | (df[COL_SLEEP] > SLEEP_MAX)
df.loc[out_of_range, COL_SLEEP] = np.nan

n_removed = n_before - df[COL_SLEEP].notna().sum()
log.info(f"  Sleep removed: {n_removed:,} ({n_removed / max(n_before, 1):.1%})")

# %% Save
cols = [COL_ID, COL_SESSION, COL_DAY, COL_START, "start_hour",
        COL_WKND, COL_HR, COL_SLEEP, COL_STEPS]
out = df[cols]
out.to_parquet(CLEAN_PARQUET, index=False)
log.info(f"Saved → {CLEAN_PARQUET}  ({len(out):,} rows)")
