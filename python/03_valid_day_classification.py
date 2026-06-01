# %% [markdown]
# # 03 — Day Validity Classification
#
# A slot is "present" (device worn) when the validity signal (default
# `min_total` = recorded minutes) exceeds `VALIDITY_MIN_VALUE`.
#
# A day is **valid** when:
# - ≥ `MIN_PRESENT_SLOTS_PER_DAY` of 12 slots are present, **and**
# - ≥ 1 present slot falls in the 10 pm – 6 am window (for sleep features).
#
# **Input:**  `data/processed/imputed/imputed_data.parquet`
# **Output:** `data/processed/day_validity.parquet`
#             `data/handoff/wear_quality.csv`

# %% Imports and setup
import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from config import (
    IMPUTED_PARQUET, DAY_VALIDITY_PARQUET, HANDOFF_DIR, LOGS_DIR,
    SESSION_LABELS,
    COL_ID, COL_SESSION, COL_DAY, COL_WKND,
    VALIDITY_SIGNAL, VALIDITY_MIN_VALUE,
    NIGHTTIME_HOURS, MIN_PRESENT_SLOTS_PER_DAY,
    MIN_VALID_DAYS, MIN_VALID_WEEKDAYS, MIN_VALID_WEEKEND_DAYS,
)

LOGS_DIR.mkdir(parents=True, exist_ok=True)
DAY_VALIDITY_PARQUET.parent.mkdir(parents=True, exist_ok=True)
HANDOFF_DIR.mkdir(parents=True, exist_ok=True)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)s  %(message)s",
    handlers=[
        logging.FileHandler(LOGS_DIR / "03_valid_day_classification.log"),
        logging.StreamHandler(),
    ],
)
log = logging.getLogger(__name__)

# %% Load imputed data
log.info(f"Loading {IMPUTED_PARQUET}")
df = pd.read_parquet(IMPUTED_PARQUET)
log.info(f"  Input rows: {len(df):,}")
assert VALIDITY_SIGNAL in df.columns, f"Validity signal '{VALIDITY_SIGNAL}' not in data"

# A slot is "present" when the wear signal exceeds the threshold
df["_present"] = df[VALIDITY_SIGNAL] > VALIDITY_MIN_VALUE

# %% Classify each (participant × session × day)
day_stats = (
    df.groupby([COL_ID, COL_SESSION, COL_DAY, COL_WKND])
    .apply(
        lambda g: pd.Series({
            "n_present_slots": int(g["_present"].sum()),
            "has_nighttime":   bool(
                g.loc[g["start_hour"].isin(NIGHTTIME_HOURS), "_present"].any()
            ),
        })
    )
    .reset_index()
)

day_stats["is_valid"] = (
    (day_stats["n_present_slots"] >= MIN_PRESENT_SLOTS_PER_DAY) &
    day_stats["has_nighttime"]
)

n_valid = day_stats["is_valid"].sum()
n_total = len(day_stats)
log.info(f"  Validity signal: {VALIDITY_SIGNAL} > {VALIDITY_MIN_VALUE}")
log.info(f"  Valid days: {n_valid:,} / {n_total:,} ({n_valid / n_total:.1%})")

# %% Aggregate per participant × session
def _summarise_wave(grp: pd.DataFrame) -> pd.Series:
    valid   = grp[grp["is_valid"]]
    n_valid = len(valid)
    n_wkday = int((~valid[COL_WKND]).sum())
    n_wkend = int(valid[COL_WKND].sum())

    valid_indices = np.where(grp["is_valid"].values)[0]
    dispersion = float(np.std(valid_indices)) if len(valid_indices) > 1 else 0.0

    n_days = len(grp)
    return pd.Series({
        "n_valid_days":         n_valid,
        "n_valid_weekdays":     n_wkday,
        "n_valid_weekend_days": n_wkend,
        "wear_time_fraction":   n_valid / n_days if n_days > 0 else 0.0,
        "valid_day_dispersion": dispersion,
    })

wave_quality = (
    day_stats.groupby([COL_ID, COL_SESSION])
    .apply(_summarise_wave)
    .reset_index()
)

# %% Apply exclusion flags
wave_quality["exclude_clustering"]       = wave_quality["n_valid_days"] < MIN_VALID_DAYS
wave_quality["exclude_weekday_features"] = wave_quality["n_valid_weekdays"] < MIN_VALID_WEEKDAYS
wave_quality["exclude_weekend_features"] = wave_quality["n_valid_weekend_days"] < MIN_VALID_WEEKEND_DAYS
wave_quality["wave_label"]               = wave_quality[COL_SESSION].map(SESSION_LABELS)

log.info(f"  Excluded from clustering:       {wave_quality['exclude_clustering'].sum()}")
log.info(f"  Excluded from weekday features: {wave_quality['exclude_weekday_features'].sum()}")
log.info(f"  Excluded from weekend features: {wave_quality['exclude_weekend_features'].sum()}")

# %% Save
day_stats.drop(columns=["_present"], errors="ignore").to_parquet(
    DAY_VALIDITY_PARQUET, index=False)
log.info(f"Saved day validity → {DAY_VALIDITY_PARQUET}")

wave_quality.to_csv(HANDOFF_DIR / "wear_quality.csv", index=False)
log.info(f"Saved wear quality → {HANDOFF_DIR / 'wear_quality.csv'}")
