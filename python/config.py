from pathlib import Path

# ---------------------------------------------------------------------------
# Active sessions
# To add the yr4 wave, append "ses-04A" to this list — no other changes needed.
# ---------------------------------------------------------------------------
SESSIONS = ["ses-02A", "ses-06A"]

SESSION_LABELS = {
    "ses-02A": "yr2",
    "ses-04A": "yr4",
    "ses-06A": "yr6",
}

# ---------------------------------------------------------------------------
# Paths (all relative to project root)
# ---------------------------------------------------------------------------
ROOT = Path(__file__).parent.parent

RAW_PARQUET        = ROOT / "dev" / "data" / "fitbit-summaries" / "activity_120m.parquet"

PROCESSED_DIR      = ROOT / "data" / "processed"
CLEAN_PARQUET      = PROCESSED_DIR / "clean"    / "clean_data.parquet"
IMPUTED_PARQUET    = PROCESSED_DIR / "imputed"  / "imputed_data.parquet"
DAY_VALIDITY_PARQUET = PROCESSED_DIR / "day_validity.parquet"
MODELS_DIR         = PROCESSED_DIR / "umap_models"

HANDOFF_DIR        = ROOT / "data" / "handoff"
LOGS_DIR           = ROOT / "logs"

# ---------------------------------------------------------------------------
# Column names (source parquet)
# ---------------------------------------------------------------------------
COL_ID      = "participant_id"
COL_SESSION = "session_id"
COL_DAY     = "day"
COL_START   = "start"
COL_WKND    = "dt_wknd"
COL_HR      = "hrate_rest_fitb"
COL_SLEEP   = "min_slp"
COL_STEPS   = "steps_total"

# ---------------------------------------------------------------------------
# Time structure
# ---------------------------------------------------------------------------
# Start hours for each 2-hour bin (e.g. hour 0 = 00:00–02:00)
ALL_SLOT_HOURS  = [0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22]
NIGHTTIME_HOURS = {22, 0, 2, 4}        # 10 pm – 6 am  (4 slots = 8 h)
DAYTIME_HOURS   = {6, 8, 10, 12, 14, 16, 18, 20}  # 6 am – 10 pm

# ---------------------------------------------------------------------------
# Artifact removal thresholds
# ---------------------------------------------------------------------------
HR_MIN          = 40    # bpm
HR_MAX          = 200   # bpm
STEPS_SPIKE_SD  = 5     # SD multiplier for within-participant-day spike detection
SLEEP_MIN       = 0     # min_slp lower bound (minutes per 2-hr slot)
SLEEP_MAX       = 120   # min_slp upper bound

# ---------------------------------------------------------------------------
# Imputation
# ---------------------------------------------------------------------------
MAX_IMPUTE_GAP       = 3   # max consecutive NaN slots to attempt imputation
MIN_NEIGHBOR_COUNT   = 3   # min valid same-hour neighbours required
NEIGHBOR_WINDOW_DAYS = 3   # ±days for same-time-of-day lookup

# ---------------------------------------------------------------------------
# Day validity
# ---------------------------------------------------------------------------
MIN_HR_SLOTS_PER_DAY    = 8   # of 12 total slots
MIN_VALID_DAYS          = 14  # of up to 21 days per wave
MIN_VALID_WEEKDAYS      = 8   # of up to 15 weekdays
MIN_VALID_WEEKEND_DAYS  = 4   # of up to 6 weekend days

# ---------------------------------------------------------------------------
# Feature extraction
# ---------------------------------------------------------------------------
MVPA_STEPS_THRESHOLD = 100   # steps per 2-hr slot to count as MVPA
L5_SLOTS  = 3                # window for L5  (~6 h at 2-hr resolution)
M10_SLOTS = 5                # window for M10 (10 h)

# ---------------------------------------------------------------------------
# Dimensionality reduction
# ---------------------------------------------------------------------------
UMAP_N_NEIGHBORS_FINAL = 30
UMAP_MIN_DIST_FINAL    = 0.1
PCA_VARIANCE_THRESHOLD = 0.95
RANDOM_STATE           = 42

# parameter grids (for sweep / documentation)
UMAP_N_NEIGHBORS = [15, 30, 50]
UMAP_MIN_DIST    = [0.0, 0.1]

# ---------------------------------------------------------------------------
# Clustering
# ---------------------------------------------------------------------------
HDBSCAN_MIN_CLUSTER_SIZES = [50, 75, 100]
HDBSCAN_MIN_SAMPLES       = [10, 20]
KMEDOIDS_K_RANGE          = range(3, 8)
N_BOOTSTRAP               = 100
BOOTSTRAP_FRACTION        = 0.8
