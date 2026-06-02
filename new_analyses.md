# Handoff Document: Adolescent Biological Rhythm Analysis
## Fitbit Cohort Study — Computational Framework Replication in R

> **Implementation:** All code lives in `project_rhythm_extraction.qmd` at the repo
> root. This document is the design/plan companion — it explains *what* each stage
> does and *why*, with the summary tables and figures the `.qmd` produces. Render
> the `.qmd` for the executable pipeline and live results.

---

## Project Overview

This project applies the computational framework from Yan & Doryab (2022) —
*"Towards a Computational Framework for Automated Discovery and Modeling of
Biological Rhythms from Wearable Data Streams"* — to a longitudinal adolescent
cohort study. The goal is to extract biological rhythm features from Fitbit wearable
data and model their development across adolescence.

**Reference paper:** Yan, R., & Doryab, A. (2022). IntelliSys 2021. LNNS 296, pp. 643–661.
https://doi.org/10.1007/978-3-030-82199-9_44

---

## Study Design

| Parameter | Detail |
|---|---|
| N participants | ~2,919 (curated cohort; the earlier "2,000" was a placeholder estimate) |
| Cohort filter | IDs with ≥2 weekends of data at **both** yr2 and yr6 (`data/ids_fitbit_v02_v06-sleep_steps_data.csv`) |
| Device | Fitbit (consumer-grade wearable) |
| Data resolution | **2-hour slots — 12 slots/day** (not minute-level as in the paper) |
| Wear duration per wave | up to ~21 days |
| Signals used | `steps_total`, `min_slp` (sleep minutes), `mets` — **HR excluded** (unreliable for this device/cohort) |
| Assessment waves | `ses-02A` (yr2 ≈ age 12), `ses-04A` (yr4 ≈ age 14), `ses-06A` (yr6 ≈ age 16) |
| Missing data | Some participants have missing days or missing slots within days |
| Raw source | `dev/data/fitbit-summaries/activity_120m.parquet` |

> **Resolution note.** Because data is aggregated to 2-hour bins, all thresholds and
> windows from the paper are re-scaled to 12-slot days (see each stage). A 24-h
> circadian period corresponds to **12 slots**; a 48-h period to **24 slots**.

---

## Analysis Pipeline

```
Raw Fitbit Data  (dev/data/fitbit-summaries/activity_120m.parquet)
      │
      ▼
[Stage 0] Data Summary  ── schema, coverage, missingness, valid-day preview
      │
      ▼
[Stage 1] Preprocessing & Imputation
      │   ├── ID filter → 3 waves
      │   ├── non-wear flag (min_total == 0)
      │   ├── valid day (≥8/12 slots) · valid wave (≥14 days)
      │   └── SMA imputation (k = ±2 slots)
      ▼
[Stage 2] Periodicity Detection  ── 4 active methods, ≥3-of-4 consensus
      │   ├── FFT
      │   ├── Chi-Squared Periodogram (Q_P)
      │   ├── Lomb-Scargle Periodogram      ← new (handles gaps natively)
      │   ├── Autocorrelation (ACF)         ← new
      │   └── Cyclic HMM approx. (depmixS4) ← DISABLED (slow); replaces Python CyHMM
      ▼
[Stage 3] Change Point Detection  ── on full wave series (~252 pts)
      │   ├── AutoNOM (BayesSpec; MAPE-selected kmax)
      │   └── Rbeast (complement + fallback)
      ▼
[Stage 4] Feature Extraction Per Participant × Wave × Signal
      │   ├── Cosinor: MESOR, amplitude, acrophase
      │   ├── Non-parametric: IS, IV, RA (M10/L5)
      │   ├── Rhythm strength: Q_P at 24 h
      │   └── Change-point: count, timing SD
      ▼
[Stage 5] Longitudinal / Growth Curve Modeling  (lme4 / lmerTest)
```

---

## Stage 0: Data Summary

Before any filtering, the `.qmd` writes summary tables straight from the raw parquet
(it never prints raw rows). These verify column names and quantify missingness so the
downstream thresholds can be sanity-checked.

**Tables produced**

| Table | Purpose |
|---|---|
| Column schema | Confirm column names/types match what the code expects |
| Participants & days per session | Coverage per wave (unfiltered) |
| Days-per-participant distribution | min / q25 / median / q75 / max days per wave |
| Signal missingness (overall + by session) | % NA and % non-wear per signal |
| Valid-day preview | % of participant-days meeting the ≥8/12-slot rule |

*Example layout (values filled at render):*

| session_id | n_participants | n_participant_days | median_days |
|---|---|---|---|
| ses-02A | … | … | … |
| ses-04A | … | … | … |
| ses-06A | … | … | … |

---

## Stage 1: Data Preprocessing

### Goals (re-scaled to 2-hr slots)
- Apply the curated ID filter and keep the 3 sessions of interest
- Non-wear = slot with `min_total == 0` (or NA) — the natural device-off indicator
- Valid **day** = ≥ 8 of 12 slots present (≤ 8 h non-wear)
- Valid **wave** = ≥ 14 valid days; otherwise drop that participant-wave
- Impute short within-wave gaps with Simple Moving Average (`imputeTS::na_ma`,
  k = ±2 slots = ±4 h); non-wear slots are set to NA first

### Flow

```
raw ──filter(ids, sessions)──► label waves ──flag non_wear──►
   day validity (≥8/12) ──► wave validity (≥14 days) ──►
   keep valid days of passing waves ──► SMA impute ──► dat_imputed
```

### Summary tables produced
- **Participant-wave counts after validity filtering** (n_total, n_pass, % pass,
  median valid days) by wave
- **Residual NA before vs. after SMA imputation** per signal

*Example layout:*

| wave | n_total | n_pass | pct_pass | median_valid_days |
|---|---|---|---|---|
| yr2 | … | … | … | … |
| yr4 | … | … | … | … |
| yr6 | … | … | … | … |

---

## Stage 2: Periodicity Detection

Run the **complementary methods** on each participant × wave × signal series, then
confirm a 24-h period if **≥ `CONSENSUS_FRAC` of the non-NA methods** detect it within
±4 h of 24 h. (The paper required all three of *its* methods; we use a proportional
threshold to stay robust at the coarser 2-h resolution.)

> **⚠️ Current config (temporary): CyHMM disabled, 4 active methods, ≥3-of-4 consensus.**
> The `depmixS4` cyclic-HMM approximation is too slow at this scale, so it is switched
> off in the `.qmd` (function + call site kept, commented). With four active methods the
> threshold is set to **`CONSENSUS_FRAC = 0.75` → ≥3 of 4 must agree** (at 0.80, 3/4 =
> 0.75 would fail and it would demand all 4). Re-enabling CyHMM restores the 5-method /
> 0.80 design below.

### The methods

| # | Method | Domain | Handles missing days | Status | Notes |
|---|---|---|---|---|---|
| 1 | FFT | Frequency | interpolated | active | dominant frequency → period |
| 2 | Chi-squared periodogram (Q_P) | variance-ratio | NA-dropped | active | Sokolove & Bushell; also gives rhythm strength |
| 3 | **Lomb-Scargle** | Frequency | **native** | active | best for gappy data; no imputation needed |
| 4 | **ACF** | lag/time | NA-aware | active | peak autocorrelation lag → period |
| 5 | **Cyclic HMM (depmixS4)** | state-transition | NA-dropped | **disabled (slow)** | BIC-selected state count → period |

> **CyHMM substitution.** The original Pierson et al. (2018) Python CyHMM
> (`github.com/epierson9/cyclic_HMMs`, file `cyclic_HMM.py`, function
> `fit_cyhmm_model()`) is **Python-2 only** and depends on a deprecated
> `pomegranate`/`scipy.misc.logsumexp` API. Rather than port it, we (a) approximate
> the cyclic HMM with `depmixS4`, and (b) add two gap-tolerant frequency methods
> (Lomb-Scargle, ACF). The pipeline stays fully in R.

### Consensus logic (conceptual)

```
per method →  period estimate (or NA)
              │
              ▼
        keep non-NA estimates                 e.g.  {FFT 24.0, Chi 24.0,
              │                                       Lomb 23.5, ACF 36.0}
              ▼
   fraction within ±4 h of 24 h  =  3/4 = 0.75  ≥ 0.75  ✔ CONFIRMED
              │                                       confirmed_period = mean(agreeing)
              ▼
   confirmed_period = mean(24.0, 24.0, 23.5) = 23.83 h
```
(With CyHMM re-enabled this becomes a 5-method / ≥4-of-5 rule at `CONSENSUS_FRAC = 0.80`.)

### Method formulas (key reference points)

- **Chi-squared Q_P** (paper Eq. 1):
  `Q_P = (K·N · Σ_h (M_h − M)²) / Σ_i (X_i − M)²`,
  significant when `Q_P > χ²(1−α, df = P−1)`, α = 0.01.
- **Acrophase / amplitude / MESOR** come later (Stage 4) once a period is confirmed.

### Summary tables & figures produced

1. **Periodicity detection results by method** — for each wave × signal, the % of
   participant-waves where each method detected 24 h, plus the final % confirmed:

   | wave | signal | n | FFT | Chi-sq | Lomb-Scargle | ACF | Confirmed (≥3/4) |
   |---|---|---|---|---|---|---|---|
   | yr2 | steps_total | … | … | … | … | … | … |
   | …  | … | … | … | … | … | … | … |

   (CyHMM column omitted while disabled.)

2. **Method-agreement bar chart** (`consensus-distribution`) — distribution of "how
   many of the 4 active methods agreed" per participant-wave, faceted by wave × signal;
   bars meeting the ≥3 threshold are highlighted. Shows *where* consensus lands and
   how often the 80 % rule is met.

3. **Confirmed-period histogram** — distribution of confirmed periods (should pile up
   at 24 h), with a dashed reference line at 24 h, faceted by wave × signal.

---

## Stage 3: Change Point Detection Within Each Period

Once a 24-h period is confirmed, detect change points (regime shifts) in the
**full wave series** (~252 points), not per day — at 2-h resolution a single day has
only 12 points, too few for the frequency-change detection AutoNOM relies on.
Restricted to participant-waves with a confirmed 24-h period.

| Method | Package | Role | Missing-data handling |
|---|---|---|---|
| **AutoNOM** | `BayesSpec` (GitHub: `hadjamar/BayesSpec`) | primary; sinusoidal RJ-MCMC segments | short gaps interpolated (≤3 days) |
| **Rbeast** | `Rbeast` (CRAN) | complement + fallback | native (`hasNA = TRUE`) |

**kmax selection (MAPE elbow, paper Table 4).** Fit AutoNOM for k ∈ {3,4,5}; pick the
k just past the largest MAPE drop. The paper found k = 4 optimal (HR MAPE 4.54 → 3.29
→ 3.14 for k = 3 → 4 → 5); re-derive on this cohort since adolescents/2-h data may
differ.

```
MAPE
  │  ●  (k=3)
  │   \
  │    ●           ← largest drop ends here ⇒ choose k = 4
  │     \____●     (k=5, marginal)
  └────────────── k
```

### Summary tables & figures produced
- **Change-point summary** — median # change points (AutoNOM & Rbeast), % with zero
  change points, by wave × signal.
- **Change-point count distribution** — dodged bar chart of CP counts per
  participant-wave, AutoNOM vs Rbeast, faceted by wave × signal.

*Example layout:*

| wave | signal | n | median_cp_autonOM | median_cp_rbeast | pct_zero_cp_autonOM | pct_zero_cp_rbeast |
|---|---|---|---|---|---|---|
| yr2 | steps_total | … | … | … | … | … |

---

## Stage 4: Feature Extraction Per Participant × Wave × Signal

For every participant × wave × signal with a confirmed 24-h period, extract:

| Family | Features | Source |
|---|---|---|
| Cosinor | MESOR, amplitude, acrophase (hour of peak) | linearized least-squares cosine fit at 24 h |
| Non-parametric circadian | IS (interdaily stability), IV (intradaily variability), M10, L5, RA (relative amplitude) | slot-level rolling windows |
| Rhythm strength | Q_P at 24 h | chi-squared periodogram |
| Change-point | n change points (AutoNOM, Rbeast), SD of CP timing | Stage 3 |
| Status | `arrhythmic` flag (no confirmed period) | Stage 2 |

**Cosinor parameterization** (period fixed at 24 h):
`y = MESOR + b·cos(2πt/24) + c·sin(2πt/24)`,
`amplitude = √(b²+c²)`, `acrophase = atan2(−c, b)` mapped to [0, 24) h.

**Non-parametric (re-scaled to slots):** M10 = highest 5-slot (10 h) mean; L5 = lowest
3-slot (6 h) mean; `RA = (M10 − L5)/(M10 + L5)`.

### Summary tables & figures produced
- **Feature summary by wave × signal** — medians of MESOR, amplitude, acrophase, IS,
  IV, RA, Q_P, plus % arrhythmic (with column spanners grouping Cosinor /
  Non-parametric / Rhythm-strength).
- **Acrophase density plot** — distribution of peak-activity hour by wave & signal
  (tests for a developmental phase delay), faceted by signal with dotted 6 h/18 h
  reference lines.

*Example layout:*

| wave | signal | n | pct_arrhythmic | med_mesor | med_amplitude | med_acrophase | med_IS | med_IV | med_RA | med_qp |
|---|---|---|---|---|---|---|---|---|---|---|
| yr2 | steps_total | … | … | … | … | … | … | … | … | … |

---

## Stage 5: Longitudinal / Growth Curve Modeling

Features at ages 12, 14, 16 → multilevel growth curves (`lme4` / `lmerTest`), one model
family per feature, fit per signal. Age centered at 12. Random intercept + slope per
participant.

```r
lmer(feature ~ age_c + (1 + age_c | participant_id), data = features_long)
# add sex / puberty when demographics are joined:
lmer(feature ~ age_c * sex + (1 + age_c | participant_id), ...)
```

### Models fit in the `.qmd`
| Outcome | Question |
|---|---|
| acrophase | Does peak-activity time shift later from 12→16? (phase delay) |
| amplitude | Does day-night contrast change with age? |
| Q_P (strength) | Does rhythm strength change with age? |
| IS, IV | Does day-to-day stability / within-day fragmentation change? |
| arrhythmic (GLMM, binomial) | Are participants more/less likely to be arrhythmic with age? |

### Key Longitudinal Research Questions
| Question | Model |
|---|---|
| Does circadian phase delay from age 12→16? | `acrophase ~ age + (age | id)` |
| Do sex differences in rhythms emerge during puberty? | `feature ~ age * sex + (age | id)` |
| Are early arrhythmic participants still arrhythmic at 16? | `arrhythmic ~ age + (1 | id)` (GLMM) |
| Does rhythm amplitude decrease with age? | `amplitude ~ age + (age | id)` |
| Does intra-day regularity (IV / SD of CP timing) change? | `IV ~ age + (age | id)` |

### Output to produce (suggested)
- Fixed-effects table per outcome (estimate, SE, t, p) — already rendered as `gt`
  tables in the `.qmd`.
- Spaghetti / predicted-trajectory plots per feature × signal (add if needed).

---

## Packages Required

```r
# CRAN packages
install.packages(c(
  "arrow", "tidyverse", "lubridate", "imputeTS",
  "lomb",        # Lomb-Scargle periodogram
  "depmixS4",    # cyclic HMM approximation
  "Rbeast",      # change-point detection
  "cosinor", "cosinor2",
  "lme4", "lmerTest", "broom.mixed",
  "knitr", "glue"   # kable tables (NOT gt — see note)
))

# GitHub packages — OPTIONAL (AutoNOM). Skipped gracefully if absent; Rbeast covers Stage 3.
# remotes::install_github("hadjamar/BayesSpec")
```

> **No `reticulate` / Python dependency** — the CyHMM step is now pure-R (`depmixS4`).
>
> **No `gt`** — tables use `knitr::kable()` instead. `gt` pulls `juicyjuice → V8`, and
> V8's wrapper will not compile under a Homebrew-GCC `~/.R/Makevars` (the
> deliberate OpenMP/`-march=native` setup): it errors with `'cmath' file not found`.
> `kable` has no compiled dependencies and matches the existing `index.qmd` convention.
>
> **`BayesSpec` is optional and renv-invisible** — referenced via `getExportedValue()`
> with a variable package name so `renv::snapshot()` does not try to install it from
> CRAN (it is GitHub-only). Install manually only if you want AutoNOM.

---

## File / Directory Structure

```
abcd_fitbit_2026/
├── dev/data/fitbit-summaries/activity_120m.parquet   # raw source (do not read directly)
├── data/
│   ├── ids_fitbit_v02_v06-sleep_steps_data.csv        # cohort ID filter
│   └── processed/features/                            # Stage 4 output (optional cache)
├── project_rhythm_extraction.qmd                      # ← the executable pipeline
├── new_analyses.md                                    # this document
└── outputs/
    ├── figures/
    └── tables/
```

---

## Known Issues and Decisions to Resolve

| Issue | Status | Resolution |
|---|---|---|
| Python CyHMM is Python-2 / deprecated deps | **Resolved** | Replaced with `depmixS4` cyclic-HMM approximation + Lomb-Scargle + ACF |
| Consensus rule (paper requires all methods) | **Resolved** | fraction-based: ≥ `CONSENSUS_FRAC` of non-NA methods within ±4 h. Design = 0.80 (≥4 of 5); **currently 0.75 (≥3 of 4) with CyHMM disabled** |
| CyHMM (depmixS4) too slow at scale | **Temporary** | Disabled for now (function + call site kept, commented); FFT + Chi-sq + Lomb-Scargle + ACF remain active |
| 2-h resolution vs minute-level paper | **Resolved** | All windows/thresholds re-scaled to 12-slot days |
| AutoNOM per-day infeasible at 12 pts/day | **Resolved** | Run AutoNOM & Rbeast on full wave series (~252 pts) |
| `BayesSpec` (AutoNOM) may not install everywhere | Open | `Rbeast` runs as complement/fallback |
| Fitbit HR noisy / unreliable | **Resolved** | HR excluded entirely |
| Non-wear vs true missing | Open | Use `min_total == 0` as non-wear proxy; refine if Fitbit wear flag available |
| kmax for AutoNOM | Open | Re-derive via MAPE elbow on this cohort (paper used k = 4) |
| Demographics (sex, puberty) for Stage 5 | Open | Join when available to enable `age * sex` models |

---

## References

- Yan, R., & Doryab, A. (2022). Towards a Computational Framework for Automated
  Discovery and Modeling of Biological Rhythms from Wearable Data Streams.
  *IntelliSys 2021*, LNNS 296, 643–661. https://doi.org/10.1007/978-3-030-82199-9_44
- Hadj-Amar, B., et al. (2019). Bayesian Model Search for Nonstationary Periodic
  Time Series. *JASA*. [AutoNOM source]
- Pierson, E., Althoff, T., & Leskovec, J. (2018). Modeling Individual Cyclic Variation
  in Human Behavior. *WWW 2018*. [CyHMM source — `github.com/epierson9/cyclic_HMMs`]
- Sokolove, P., & Bushell, W. (1978). The chi square periodogram. *J. Theor. Biol.*
- Lomb, N. R. (1976); Scargle, J. D. (1982). Least-squares spectral analysis for
  unevenly-spaced data. [Lomb-Scargle source]
- Visser, I., & Speekenbrink, M. (2010). depmixS4: An R Package for Hidden Markov
  Models. *J. Stat. Softw.*
- Zhao, K., et al. (2019). Rbeast: Bayesian change-point detection and time series
  decomposition.
- Moritz, S., & Bartz-Beielstein, T. (2017). imputeTS. *R Journal*.
```
