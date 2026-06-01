# %% [markdown]
# # 02 — Within-Day Imputation
#
# Imputes gaps of ≤ 3 consecutive missing 2-hr slots using the
# same-time-of-day median from ± 3 calendar days (≥ 3 neighbours required).
# Runs on every signal in `config.SIGNALS`. The validity signal (e.g. min_total)
# is carried through untouched — it is a wear indicator, not imputed.
#
# **Input:**  `data/processed/clean/clean_data.parquet`
# **Output:** `data/processed/imputed/imputed_data.parquet`

# %% Imports and setup
import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from config import (
    CLEAN_PARQUET, IMPUTED_PARQUET, LOGS_DIR,
    SIGNAL_NAMES, VALIDITY_SIGNAL,
    COL_ID, COL_SESSION, COL_DAY, COL_WKND,
    ALL_SLOT_HOURS,
    MAX_IMPUTE_GAP, MIN_NEIGHBOR_COUNT, NEIGHBOR_WINDOW_DAYS,
)

LOGS_DIR.mkdir(parents=True, exist_ok=True)
IMPUTED_PARQUET.parent.mkdir(parents=True, exist_ok=True)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)s  %(message)s",
    handlers=[
        logging.FileHandler(LOGS_DIR / "02_imputation.log"),
        logging.StreamHandler(),
    ],
)
log = logging.getLogger(__name__)

# Columns carried through but not imputed (wear / validity indicators)
PASSTHROUGH = [VALIDITY_SIGNAL] if VALIDITY_SIGNAL not in SIGNAL_NAMES else []

# %% Load clean data
log.info(f"Loading {CLEAN_PARQUET}")
df = pd.read_parquet(CLEAN_PARQUET)
log.info(f"  Input rows: {len(df):,}  signals={SIGNAL_NAMES}")

for col in SIGNAL_NAMES:
    assert col in df.columns, f"Missing column: {col}"

df["start_hour"] = df["start_hour"].astype(int)

# %% Build complete time grid
# For every (participant, session, day) observed, ensure all 12 slots exist.
# Slots absent from the source data are treated as NaN (and remain non-present).

day_index = df[[COL_ID, COL_SESSION, COL_DAY, COL_WKND]].drop_duplicates()

grid = pd.DataFrame(
    [
        {COL_ID: r[COL_ID], COL_SESSION: r[COL_SESSION],
         COL_DAY: r[COL_DAY], COL_WKND: r[COL_WKND], "start_hour": h}
        for _, r in day_index.iterrows()
        for h in ALL_SLOT_HOURS
    ]
)

value_cols = SIGNAL_NAMES + PASSTHROUGH
merged = grid.merge(
    df[[COL_ID, COL_SESSION, COL_DAY, "start_hour"] + value_cols],
    on=[COL_ID, COL_SESSION, COL_DAY, "start_hour"],
    how="left",
)
log.info(f"  Complete grid: {len(merged):,} rows "
         f"({len(merged) - len(df):,} slots added as NaN)")

# %% Imputation helpers

def _imputable_mask(nan_mask: np.ndarray, max_gap: int) -> np.ndarray:
    """Boolean mask of NaN positions that fall within a run of length ≤ max_gap."""
    n = len(nan_mask)
    out = np.zeros(n, dtype=bool)
    i = 0
    while i < n:
        if nan_mask[i]:
            j = i
            while j < n and nan_mask[j]:
                j += 1
            if j - i <= max_gap:
                out[i:j] = True
            i = j
        else:
            i += 1
    return out


def _impute_group(group: pd.DataFrame) -> pd.DataFrame:
    """Impute all configured signals for one (participant, session) group."""
    group = group.sort_values([COL_DAY, "start_hour"]).copy()

    sorted_days = sorted(group[COL_DAY].unique())
    n_days = len(sorted_days)

    for signal in SIGNAL_NAMES:
        pivot = (
            group.pivot(index=COL_DAY, columns="start_hour", values=signal)
            .reindex(index=sorted_days, columns=ALL_SLOT_HOURS)
        )
        arr = pivot.values.astype(float)   # (n_days × 12)

        for di in range(n_days):
            row = arr[di]
            nan_mask = np.isnan(row)
            if not nan_mask.any():
                continue

            imputable = _imputable_mask(nan_mask, MAX_IMPUTE_GAP)
            for si in np.where(imputable)[0]:
                neighbours = [
                    arr[di + delta, si]
                    for delta in range(-NEIGHBOR_WINDOW_DAYS, NEIGHBOR_WINDOW_DAYS + 1)
                    if delta != 0
                    and 0 <= di + delta < n_days
                    and not np.isnan(arr[di + delta, si])
                ]
                if len(neighbours) >= MIN_NEIGHBOR_COUNT:
                    arr[di, si] = float(np.median(neighbours))

        # Write imputed values back (vectorized via MultiIndex)
        idx = pd.MultiIndex.from_product(
            [sorted_days, ALL_SLOT_HOURS], names=[COL_DAY, "start_hour"]
        )
        group = group.set_index([COL_DAY, "start_hour"])
        group[signal] = pd.Series(arr.flatten(), index=idx)
        group = group.reset_index()

    return group

# %% Run imputation
log.info("Running imputation...")

parts = []
groups = list(merged.groupby([COL_ID, COL_SESSION]))
n_groups = len(groups)

for i, ((pid, sess), grp) in enumerate(groups):
    if i % 500 == 0:
        log.info(f"  {i + 1}/{n_groups}")
    parts.append(_impute_group(grp))

imputed = pd.concat(parts, ignore_index=True)

# %% Report imputed counts
for signal in SIGNAL_NAMES:
    n_imputed = imputed[signal].notna().sum() - merged[signal].notna().sum()
    log.info(f"  {signal}: {n_imputed:,} slots imputed")

# %% Save
imputed.to_parquet(IMPUTED_PARQUET, index=False)
log.info(f"Saved → {IMPUTED_PARQUET}  ({len(imputed):,} rows)")
