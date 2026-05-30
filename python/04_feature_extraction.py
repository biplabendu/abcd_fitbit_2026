# %% [markdown]
# # 04 — Feature Extraction
#
# Computes circadian, HR, activity, and sleep features from valid days only.
# Features are z-scored at the population level within each session.
#
# **Input:**  `data/processed/imputed/imputed_data.parquet`
#             `data/processed/day_validity.parquet`
#             `data/handoff/wear_quality.csv`
# **Output:** `data/handoff/features_yr2.csv`
#             `data/handoff/features_yr6.csv`

# %% Imports and setup
import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from config import (
    IMPUTED_PARQUET, DAY_VALIDITY_PARQUET, HANDOFF_DIR, LOGS_DIR,
    SESSIONS, SESSION_LABELS,
    COL_ID, COL_SESSION, COL_DAY, COL_WKND,
    COL_HR, COL_SLEEP, COL_STEPS,
    ALL_SLOT_HOURS, NIGHTTIME_HOURS, DAYTIME_HOURS,
    MVPA_STEPS_THRESHOLD, L5_SLOTS, M10_SLOTS,
)
from utils.circadian import fit_cosinor, interdaily_stability, intradaily_variability, l5_m10

LOGS_DIR.mkdir(parents=True, exist_ok=True)
HANDOFF_DIR.mkdir(parents=True, exist_ok=True)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)s  %(message)s",
    handlers=[
        logging.FileHandler(LOGS_DIR / "04_feature_extraction.log"),
        logging.StreamHandler(),
    ],
)
log = logging.getLogger(__name__)

SLOT_HOURS_ARR = np.array(ALL_SLOT_HOURS, dtype=float)

# %% Load and merge data
log.info("Loading data...")
df          = pd.read_parquet(IMPUTED_PARQUET)
day_valid   = pd.read_parquet(DAY_VALIDITY_PARQUET)
wear        = pd.read_csv(HANDOFF_DIR / "wear_quality.csv")

df = df.merge(
    day_valid[[COL_ID, COL_SESSION, COL_DAY, "is_valid"]],
    on=[COL_ID, COL_SESSION, COL_DAY], how="left",
)
df = df.merge(
    wear[[COL_ID, COL_SESSION,
          "exclude_weekend_features", "exclude_weekday_features", "exclude_clustering"]],
    on=[COL_ID, COL_SESSION], how="left",
)

valid = df[df["is_valid"] == True].copy()
log.info(f"  Valid-day rows: {len(valid):,}")

# %% Helpers

def _autocorr_lag1(series: pd.Series) -> float:
    s = series.dropna()
    return float(s.autocorr(lag=1)) if len(s) >= 4 else np.nan


def _circ_mean(hours: pd.Series) -> float:
    v = hours.dropna()
    if len(v) == 0:
        return np.nan
    r = v * (2 * np.pi / 24)
    return float(np.arctan2(np.sin(r).mean(), np.cos(r).mean()) * 24 / (2 * np.pi) % 24)


def _circ_sd(hours: pd.Series) -> float:
    v = hours.dropna()
    if len(v) < 2:
        return np.nan
    r = v * (2 * np.pi / 24)
    R = np.sqrt(np.sin(r).mean() ** 2 + np.cos(r).mean() ** 2)
    return float(np.sqrt(-2 * np.log(np.clip(R, 1e-9, 1))) * 24 / (2 * np.pi))

# %% Extract features per participant × session
all_rows = []
groups   = list(valid.groupby([COL_ID, COL_SESSION]))
n_groups = len(groups)

for i, ((pid, sess), grp) in enumerate(groups):
    if i % 500 == 0:
        log.info(f"  {i + 1}/{n_groups}")

    row = {COL_ID: pid, COL_SESSION: sess,
           "wave_label": SESSION_LABELS[sess],
           "n_valid_days": grp[COL_DAY].nunique()}

    excl_wknd  = bool(grp["exclude_weekend_features"].iloc[0])
    excl_wkday = bool(grp["exclude_weekday_features"].iloc[0])

    # --------------------------------------------------------
    # 2.1  Circadian features
    # --------------------------------------------------------
    hr_profile    = grp.groupby("start_hour")[COL_HR].mean().reindex(ALL_SLOT_HOURS).values
    steps_profile = grp.groupby("start_hour")[COL_STEPS].mean().reindex(ALL_SLOT_HOURS).values

    _, row["hr_amplitude"],       row["hr_acrophase"]       = fit_cosinor(SLOT_HOURS_ARR, hr_profile)
    _, row["activity_amplitude"], row["activity_acrophase"] = fit_cosinor(SLOT_HOURS_ARR, steps_profile)

    steps_2d = (
        grp.pivot_table(index=COL_DAY, columns="start_hour",
                        values=COL_STEPS, aggfunc="mean")
        .reindex(columns=ALL_SLOT_HOURS).values.astype(float)
    )
    row["interdaily_stability"]   = interdaily_stability(steps_2d)
    row["intradaily_variability"] = intradaily_variability(steps_2d)

    lm = l5_m10(steps_profile, ALL_SLOT_HOURS, L5_SLOTS, M10_SLOTS)
    row.update({"L5": lm["l5"], "L5_onset": lm["l5_onset"],
                "M10": lm["m10"], "M10_onset": lm["m10_onset"],
                "relative_amplitude": lm["relative_amplitude"]})

    # --------------------------------------------------------
    # 2.2  HR summary
    # --------------------------------------------------------
    daily_hr          = grp.groupby(COL_DAY)[COL_HR].mean()
    row["hr_mean_daily"]    = float(daily_hr.mean())
    row["hr_sd_daily"]      = float(daily_hr.std())
    row["hr_cv"]            = row["hr_sd_daily"] / row["hr_mean_daily"] if row["hr_mean_daily"] > 0 else np.nan
    row["hr_autocorr_lag1"] = _autocorr_lag1(daily_hr)

    daytime = grp[grp["start_hour"].isin(DAYTIME_HOURS)]
    row["hr_active_mean"]  = float(daytime.loc[daytime[COL_STEPS] > 0,  COL_HR].mean())
    row["hr_resting_mean"] = float(daytime.loc[daytime[COL_STEPS] == 0, COL_HR].mean())

    # --------------------------------------------------------
    # 2.3  Steps summary
    # --------------------------------------------------------
    daily_steps            = grp.groupby(COL_DAY)[COL_STEPS].mean()
    row["steps_mean_daily"]    = float(daily_steps.mean())
    row["steps_sd_daily"]      = float(daily_steps.std())
    row["steps_cv"]            = row["steps_sd_daily"] / row["steps_mean_daily"] if row["steps_mean_daily"] > 0 else np.nan
    row["steps_autocorr_lag1"] = _autocorr_lag1(daily_steps)

    dt_slots = daytime[COL_STEPS]
    n_dt     = dt_slots.notna().sum()
    row["sedentary_fraction"] = float((dt_slots == 0).sum() / n_dt) if n_dt > 0 else np.nan
    row["mvpa_fraction"]      = float((dt_slots >= MVPA_STEPS_THRESHOLD).sum() / n_dt) if n_dt > 0 else np.nan

    if not excl_wknd:
        wkday_steps = grp.loc[~grp[COL_WKND], COL_STEPS].mean()
        wkend_steps = grp.loc[grp[COL_WKND],  COL_STEPS].mean()
        row["weekend_weekday_delta"] = float(wkend_steps - wkday_steps)
    else:
        row["weekend_weekday_delta"] = np.nan

    # --------------------------------------------------------
    # 2.4  Sleep features (approximated from min_slp)
    # --------------------------------------------------------
    night = grp[grp["start_hour"].isin(NIGHTTIME_HOURS)]

    nightly_dur = night.groupby(COL_DAY)[COL_SLEEP].sum() / 60.0   # hours
    row["sleep_duration_mean"]  = float(nightly_dur.mean())
    row["sleep_duration_sd"]    = float(nightly_dur.std())
    row["sleep_autocorr_lag1"]  = _autocorr_lag1(nightly_dur)

    def _onset(day_grp):
        s = day_grp.loc[day_grp[COL_SLEEP] > 0].sort_values("start_hour")
        return float(s["start_hour"].iloc[0]) if len(s) > 0 else np.nan

    nightly_onset = night.groupby(COL_DAY).apply(_onset)
    row["sleep_onset_mean"] = _circ_mean(nightly_onset)
    row["sleep_onset_sd"]   = _circ_sd(nightly_onset)

    night_sleeping         = night[night[COL_SLEEP] > 0]
    row["sleep_efficiency_mean"] = float((night_sleeping[COL_SLEEP] / 120).mean()) if len(night_sleeping) > 0 else np.nan

    if not excl_wknd:
        def _midpoint(day_grp):
            s = day_grp.loc[day_grp[COL_SLEEP] > 0].sort_values("start_hour")
            if len(s) == 0:
                return np.nan
            return (s["start_hour"].iloc[0] + s["start_hour"].iloc[-1] + 2) / 2.0

        nightly_mid  = night.groupby(COL_DAY).apply(_midpoint)
        wkday_flag   = grp.groupby(COL_DAY)[COL_WKND].first()
        wkday_mid    = nightly_mid[~wkday_flag].mean()
        wkend_mid    = nightly_mid[wkday_flag].mean()
        row["social_jet_lag"] = float(wkend_mid - wkday_mid)
    else:
        row["social_jet_lag"] = np.nan

    all_rows.append(row)

# %% Assemble feature matrix
features = pd.DataFrame(all_rows)
feature_cols = [
    c for c in features.columns
    if c not in [COL_ID, COL_SESSION, "wave_label", "n_valid_days"]
]
log.info(f"  Feature columns: {len(feature_cols)}")

# %% Population-level z-score within each session
for sess in SESSIONS:
    mask = features[COL_SESSION] == sess
    for col in feature_cols:
        mu = features.loc[mask, col].mean()
        sd = features.loc[mask, col].std()
        features.loc[mask, f"z_{col}"] = (
            (features.loc[mask, col] - mu) / sd if sd > 0 else 0.0
        )

# %% Save one CSV per wave
for sess in SESSIONS:
    label   = SESSION_LABELS[sess]
    wave_df = features[features[COL_SESSION] == sess]
    path    = HANDOFF_DIR / f"features_{label}.csv"
    wave_df.to_csv(path, index=False)
    log.info(f"  Saved {label}: {len(wave_df):,} participants → {path}")

log.info("Feature extraction complete.")
