# %% [markdown]
# # 04 — Feature Extraction
#
# Computes features from valid days only, per participant per wave.
# Feature families are driven by each signal's `kind` in `config.SIGNALS`:
#
# - **activity** signals → cosinor (amplitude, acrophase), interdaily stability,
#   intradaily variability, L5/M10, daily mean/sd/cv/autocorr, and (if a
#   `mvpa_threshold` is set) sedentary / MVPA fractions and weekend–weekday delta.
# - **sleep** signals → nightly duration mean/sd, onset mean/sd (circular),
#   efficiency, social jet lag, and day-to-day autocorrelation.
#
# All feature columns are prefixed by their signal name (e.g. `steps_total_M10`).
# Features are population z-scored within each session.
#
# **Input:**  `data/processed/imputed/imputed_data.parquet`,
#             `data/processed/day_validity.parquet`,
#             `data/handoff/wear_quality.csv`
# **Output:** `data/handoff/features_{label}.csv` (one per session)
#             `data/handoff/diurnal_profiles.csv` (per-participant mean diurnal
#             profile per signal × slot, for the R diurnal-curve figure)

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
    SIGNALS, SIGNAL_NAMES, ACTIVITY_SIGNALS, SLEEP_SIGNALS,
    COL_ID, COL_SESSION, COL_DAY, COL_WKND,
    ALL_SLOT_HOURS, NIGHTTIME_HOURS, DAYTIME_HOURS,
    L5_SLOTS, M10_SLOTS,
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

SLOT_HOURS_ARR   = np.array(ALL_SLOT_HOURS, dtype=float)
# Sort order within a midnight-wrapping night: evening (22) first, then early morning
NIGHT_SLOT_ORDER = {22: 0, 0: 1, 2: 2, 4: 3}

# %% Load and merge data
log.info("Loading data...")
df        = pd.read_parquet(IMPUTED_PARQUET)
day_valid = pd.read_parquet(DAY_VALIDITY_PARQUET)
wear      = pd.read_csv(HANDOFF_DIR / "wear_quality.csv")

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
log.info(f"  Valid-day rows: {len(valid):,}  "
         f"(activity={ACTIVITY_SIGNALS}, sleep={SLEEP_SIGNALS})")

# %% Helpers

def _autocorr_lag1(series: pd.Series) -> float:
    s = series.dropna()
    if len(s) < 4 or s.std() == 0:
        return np.nan
    return float(s.autocorr(lag=1))


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


def _annotate_night_anchor(night_df: pd.DataFrame) -> pd.DataFrame:
    """Tag each nighttime row with a midnight-wrapping night anchor.

    Hour-22 slots anchor to their own calendar day value.  Hours 0/2/4 are
    the early-morning tail of the *previous* night, so they are reassigned
    to the immediately preceding day value in the sorted day sequence.
    Works regardless of whether COL_DAY is a date string, Timestamp, or
    integer index — only sortability is required.

    Also adds night_slot_order (22→0, 0→1, 2→2, 4→3) for correct onset
    sorting across the midnight boundary.
    """
    sorted_days = sorted(night_df[COL_DAY].unique())
    prev_day    = {d: sorted_days[i - 1] if i > 0 else d
                   for i, d in enumerate(sorted_days)}

    is_early = night_df["start_hour"].isin({0, 2, 4})
    anchors  = night_df[COL_DAY].copy()
    anchors.loc[is_early] = night_df.loc[is_early, COL_DAY].map(prev_day).values

    return night_df.assign(
        night_anchor     = anchors,
        night_slot_order = night_df["start_hour"].map(NIGHT_SLOT_ORDER),
    )


def _activity_features(grp: pd.DataFrame, sig: str, meta: dict, excl_wknd: bool) -> dict:
    """Circadian + daily-summary features for one activity-type signal."""
    out = {}

    # Mean diurnal profile → cosinor + L5/M10
    profile = grp.groupby("start_hour")[sig].mean().reindex(ALL_SLOT_HOURS).values
    _, out[f"{sig}_amplitude"], out[f"{sig}_acrophase"] = fit_cosinor(SLOT_HOURS_ARR, profile)

    mat = (
        grp.pivot_table(index=COL_DAY, columns="start_hour", values=sig, aggfunc="mean")
        .reindex(columns=ALL_SLOT_HOURS).values.astype(float)
    )
    out[f"{sig}_IS"] = interdaily_stability(mat)
    out[f"{sig}_IV"] = intradaily_variability(mat)

    lm = l5_m10(profile, ALL_SLOT_HOURS, L5_SLOTS, M10_SLOTS)
    out[f"{sig}_L5"]                 = lm["l5"]
    out[f"{sig}_L5_onset"]           = lm["l5_onset"]
    out[f"{sig}_M10"]                = lm["m10"]
    out[f"{sig}_M10_onset"]          = lm["m10_onset"]
    out[f"{sig}_relative_amplitude"] = lm["relative_amplitude"]

    # Daily summaries
    daily = grp.groupby(COL_DAY)[sig].mean()
    mean_ = float(daily.mean())
    sd_   = float(daily.std())
    out[f"{sig}_mean_daily"]    = mean_
    out[f"{sig}_sd_daily"]      = sd_
    out[f"{sig}_cv"]            = sd_ / mean_ if mean_ > 0 else np.nan
    out[f"{sig}_autocorr_lag1"] = _autocorr_lag1(daily)

    # Sedentary / MVPA fractions (daytime slots only)
    thr = meta.get("mvpa_threshold")
    if thr is not None:
        dt = grp.loc[grp["start_hour"].isin(DAYTIME_HOURS), sig]
        n_dt = dt.notna().sum()
        out[f"{sig}_sedentary_fraction"] = float((dt == 0).sum() / n_dt) if n_dt > 0 else np.nan
        out[f"{sig}_mvpa_fraction"]      = float((dt >= thr).sum() / n_dt) if n_dt > 0 else np.nan

    # Weekend–weekday delta
    if not excl_wknd:
        wk = grp.loc[~grp[COL_WKND], sig].mean()
        we = grp.loc[grp[COL_WKND],  sig].mean()
        out[f"{sig}_weekend_weekday_delta"] = float(we - wk)
    else:
        out[f"{sig}_weekend_weekday_delta"] = np.nan

    return out


def _sleep_features(grp: pd.DataFrame, sig: str, slot_minutes: float, excl_wknd: bool) -> dict:
    """Nightly duration / timing / efficiency features for one sleep-type signal.

    Nights are anchored to the hour-22 slot (evening) and extend forward into
    the following day's 0/2/4 slots, so each biological night is grouped
    correctly across the midnight boundary.
    """
    out = {}
    night = _annotate_night_anchor(grp[grp["start_hour"].isin(NIGHTTIME_HOURS)].copy())

    nightly_dur = night.groupby("night_anchor")[sig].sum() / 60.0   # hours
    out[f"{sig}_duration_mean"] = float(nightly_dur.mean())
    out[f"{sig}_duration_sd"]   = float(nightly_dur.std())
    out[f"{sig}_autocorr_lag1"] = _autocorr_lag1(nightly_dur)

    def _onset(night_grp):
        s = night_grp.loc[night_grp[sig] > 0].sort_values("night_slot_order")
        return float(s["start_hour"].iloc[0]) if len(s) > 0 else np.nan

    nightly_onset = night.groupby("night_anchor").apply(_onset)
    out[f"{sig}_onset_mean"] = _circ_mean(nightly_onset)
    out[f"{sig}_onset_sd"]   = _circ_sd(nightly_onset)

    sleeping = night[night[sig] > 0]
    out[f"{sig}_efficiency_mean"] = float((sleeping[sig] / slot_minutes).mean()) if len(sleeping) > 0 else np.nan

    if not excl_wknd:
        def _midpoint(night_grp):
            s = night_grp.loc[night_grp[sig] > 0].sort_values("night_slot_order")
            if len(s) == 0:
                return np.nan
            return (s["start_hour"].iloc[0] + s["start_hour"].iloc[-1] + 2) / 2.0

        mid = night.groupby("night_anchor").apply(_midpoint)

        # SJL = free-night midsleep − work-night midsleep (conventional sign: positive
        # when sleeping later on free nights).  Weekend status comes from the hour-22
        # (evening) slot that anchors each night.
        evening     = night[night["start_hour"] == 22]
        is_free_night = evening.groupby("night_anchor")[COL_WKND].first().reindex(mid.index)
        out[f"{sig}_social_jet_lag"] = float(
            mid[is_free_night == True].mean() - mid[is_free_night == False].mean()
        )
    else:
        out[f"{sig}_social_jet_lag"] = np.nan

    return out

# %% Extract features per participant × session
all_rows     = []
diurnal_rows = []   # per-participant mean diurnal profile per signal × slot
groups   = list(valid.groupby([COL_ID, COL_SESSION]))
n_groups = len(groups)

for i, ((pid, sess), grp) in enumerate(groups):
    if i % 500 == 0:
        log.info(f"  {i + 1}/{n_groups}")

    row = {COL_ID: pid, COL_SESSION: sess,
           "wave_label": SESSION_LABELS[sess],
           "n_valid_days": grp[COL_DAY].nunique()}

    excl_wknd = bool(grp["exclude_weekend_features"].iloc[0])

    for sig in ACTIVITY_SIGNALS:
        row.update(_activity_features(grp, sig, SIGNALS[sig], excl_wknd))

    for sig in SLEEP_SIGNALS:
        slot_minutes = SIGNALS[sig].get("clip_max", 120) or 120
        row.update(_sleep_features(grp, sig, slot_minutes, excl_wknd))

    # Mean diurnal profile (mean across valid days per slot) for every signal
    for sig in SIGNAL_NAMES:
        profile = grp.groupby("start_hour")[sig].mean()
        for hour in ALL_SLOT_HOURS:
            diurnal_rows.append({
                COL_ID: pid, COL_SESSION: sess,
                "wave_label": SESSION_LABELS[sess],
                "signal": sig, "start_hour": hour,
                "value": float(profile.get(hour, np.nan)),
            })

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

# %% Save diurnal profiles (for the R diurnal-curve figure, Fig 2)
diurnal = pd.DataFrame(diurnal_rows)
diurnal_path = HANDOFF_DIR / "diurnal_profiles.csv"
diurnal.to_csv(diurnal_path, index=False)
log.info(f"  Saved diurnal profiles: {len(diurnal):,} rows → {diurnal_path}")

log.info("Feature extraction complete.")
