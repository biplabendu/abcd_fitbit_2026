# Handoff Document: Adolescent Biological Rhythm Analysis
## Fitbit Cohort Study — Computational Framework Replication in R

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
| N participants | 2,000 |
| Device | Fitbit (consumer-grade wearable) |
| Wear duration per wave | 20 days |
| Signals available | Heart rate (1/min), steps, sleep stages; possibly skin temperature |
| Assessment waves | Ages 12, 14, and 16 (3 waves total) |
| Missing data | Some participants have missing days or missing hours within days |

---

## Analysis Pipeline

The pipeline follows three sequential stages, matching the paper:

```
Raw Fitbit Data
      │
      ▼
[Stage 1] Data Preprocessing & Imputation
      │
      ▼
[Stage 2] Periodicity Detection
      │   ├── Fast Fourier Transform (FFT)
      │   ├── Chi-Squared Periodogram
      │   └── Cyclic Hidden Markov Models (CyHMMs)
      │
      ▼
[Stage 3] Change Point Detection Within Each Period
      │   └── AutoNOM (or Rbeast as fallback)
      │
      ▼
[Stage 4] Feature Extraction Per Participant Per Wave
      │
      ▼
[Stage 5] Longitudinal / Growth Curve Modeling
```

---

## Stage 1: Data Preprocessing

### Goals
- Load per-participant Fitbit data (one file or row-set per participant per wave)
- Distinguish non-wear from true missing (non-wear = consecutive zeros or gaps > 60 min)
- Impute missing hours within days using Simple Moving Average (SMA)
- Flag and optionally exclude days with > 4 hours of non-wear
- Exclude participants with fewer than 14 valid days in a wave

### Expected Input Format

```
participant_id | wave | datetime            | heart_rate | steps | sleep_stage
P001           | 1    | 2023-01-01 00:00:00 | 58         | 0     | 3
P001           | 1    | 2023-01-01 00:01:00 | 59         | 0     | 3
...
```

### R Implementation

```r
library(tidyverse)
library(imputeTS)
library(lubridate)

preprocess_participant <- function(df, min_valid_days = 14, max_missing_hours = 4) {

  df <- df %>%
    arrange(datetime) %>%
    mutate(
      date = as_date(datetime),
      hour = hour(datetime)
    )

  # Flag non-wear: gaps > 60 consecutive minutes with zero HR
  df <- df %>%
    mutate(
      non_wear = (heart_rate == 0 | is.na(heart_rate))
    )

  # Count valid hours per day
  valid_days <- df %>%
    group_by(date) %>%
    summarise(missing_hours = sum(non_wear) / 60) %>%
    filter(missing_hours <= max_missing_hours) %>%
    pull(date)

  if (length(valid_days) < min_valid_days) {
    return(NULL)  # exclude participant-wave
  }

  df <- df %>% filter(date %in% valid_days)

  # SMA imputation for remaining missing minutes
  df <- df %>%
    mutate(
      heart_rate = na_ma(heart_rate, k = 30, weighting = "simple"),
      steps      = na_ma(steps,      k = 30, weighting = "simple")
    )

  return(df)
}
```

---

## Stage 2: Periodicity Detection

Run all three methods on each participant-wave-signal combination.
A period is considered confirmed if detected by at least two of the three methods
at p < 0.01.

### 2a. Fast Fourier Transform (FFT)

```r
run_fft <- function(x, sampling_interval_min = 1) {
  n <- length(x)
  fft_result <- fft(x)
  power <- Mod(fft_result)^2
  freq <- (0:(n - 1)) / n  # cycles per sample

  # Convert to cycles per hour
  freq_per_hour <- freq / (sampling_interval_min / 60)

  # Exclude DC component (index 1) and upper half (mirror)
  usable <- 2:(n %/% 2)
  dominant_idx <- usable[which.max(power[usable])]
  dominant_period_hours <- 1 / freq_per_hour[dominant_idx]

  return(list(
    dominant_period = dominant_period_hours,
    power_spectrum  = data.frame(freq = freq_per_hour[usable],
                                 power = power[usable])
  ))
}
```

### 2b. Chi-Squared Periodogram

Implements the Sokolove & Bushell Q_P statistic from the paper (Eq. 1).

```r
chi_sq_periodogram <- function(x, period_range = c(12, 48),
                                step = 0.1, alpha = 0.01) {
  N <- length(x)
  M <- mean(x, na.rm = TRUE)
  total_var <- sum((x - M)^2, na.rm = TRUE)
  periods <- seq(period_range[1], period_range[2], by = step)

  qp_values <- sapply(periods, function(P) {
    P_int <- round(P)
    K <- floor(N / P_int)
    if (K < 2) return(NA)
    segments <- matrix(x[1:(K * P_int)], nrow = P_int, ncol = K)
    Mh <- rowMeans(segments, na.rm = TRUE)
    QP <- (K * N * sum((Mh - M)^2)) / total_var
    return(QP)
  })

  # Significance threshold: chi-sq with P-1 df
  sig_threshold <- qchisq(1 - alpha, df = round(periods) - 1)
  significant <- qp_values > sig_threshold

  return(data.frame(
    period   = periods,
    QP       = qp_values,
    sig_threshold = sig_threshold,
    significant = significant
  ))
}
```

### 2c. Cyclic Hidden Markov Models (CyHMMs)

CyHMMs require constrained cyclic transition matrices. Use `depmixS4` with
manual constraints, or call the original Python implementation via `reticulate`.

**Option A — depmixS4 (approximate, constrained HMM):**

```r
library(depmixS4)

run_cyhmm <- function(x, n_states_range = 20:30) {
  # Try different cycle lengths (n_states = candidate period in hours
  # assuming hourly aggregated data)
  results <- lapply(n_states_range, function(n_states) {
    tryCatch({
      mod <- depmix(response = list(x ~ 1),
                    data = data.frame(x = x),
                    nstates = n_states,
                    family = list(gaussian()))
      # Fix transition matrix to be cyclic (each state -> next state only)
      # This requires custom constraint setup — see depmixS4 vignette
      fit <- fit(mod, verbose = FALSE)
      return(list(period = n_states, BIC = BIC(fit)))
    }, error = function(e) NULL)
  })

  best <- results[[which.min(sapply(results, function(r) r$BIC))]]
  return(best$period)
}
```

**Option B — Python CyHMM via reticulate (preferred for fidelity):**

```r
library(reticulate)
# use_python("/usr/bin/python3")
# source_python("cyhmm_pierson2018.py")  # obtain from Pierson et al. 2018 repo
# period <- run_cyhmm_python(x)
```

### Period Confirmation Logic

```r
confirm_period <- function(fft_period, chi_periods, cyhmm_period,
                           tolerance_hours = 2) {
  candidates <- chi_periods$period[chi_periods$significant]
  fft_match  <- any(abs(candidates - fft_period) < tolerance_hours)
  cyhmm_match <- any(abs(candidates - cyhmm_period) < tolerance_hours)

  if (fft_match & cyhmm_match) return(fft_period)
  if (fft_match | cyhmm_match) return(fft_period)
  return(NA)  # no consensus
}
```

---

## Stage 3: Change Point Detection Within Each Period

Once a 24-hour period is confirmed, apply AutoNOM (or Rbeast as fallback)
to each day's data to find intra-day change points.

### 3a. AutoNOM (Primary — requires BayesSpec from GitHub)

```r
# devtools::install_github("hadjamar/BayesSpec")
library(BayesSpec)

run_autonon_day <- function(day_vector, kmax = 4, mmax = 3) {
  result <- autnom(y = day_vector, kmax = kmax, mmax = mmax)
  return(list(
    change_points = result$changepoints,
    n_segments    = length(result$changepoints) + 1,
    fitted        = result$fitted
  ))
}

# Select optimal kmax using MAPE
select_kmax <- function(day_vector, k_range = 3:5) {
  mape_vals <- sapply(k_range, function(k) {
    res <- run_autonon_day(day_vector, kmax = k)
    actual    <- day_vector
    predicted <- res$fitted
    mean(abs((predicted - actual) / actual), na.rm = TRUE) * 100
  })
  optimal_k <- k_range[which.min(diff(mape_vals)) + 1]  # elbow
  return(optimal_k)
}
```

### 3b. Rbeast (Fallback — available on CRAN)

```r
library(Rbeast)

run_beast_day <- function(day_vector, period_in_samples = 1440) {
  result <- beast(
    day_vector,
    season  = "harmonic",
    period  = period_in_samples,
    quiet   = TRUE
  )
  change_points <- which(result$trend$cpPr > 0.5)
  return(list(
    change_points = change_points,
    fitted        = result$trend$Y
  ))
}
```

---

## Stage 4: Feature Extraction Per Participant Per Wave

After running Stages 2 and 3, extract a feature vector for each participant × wave.

```r
extract_rhythm_features <- function(participant_id, wave, signal_vec,
                                    confirmed_period, cp_results) {

  # Cosinor fit for MESOR, amplitude, acrophase
  library(cosinor)
  t <- seq_along(signal_vec) / 60  # time in hours
  cos_fit <- cosinor.lm(signal_vec ~ time(t), period = confirmed_period,
                         data = data.frame(signal_vec, t))
  cosinor_params <- summary(cos_fit)

  # Change point features across all days
  all_cps <- unlist(lapply(cp_results, `[[`, "change_points"))

  data.frame(
    participant_id        = participant_id,
    wave                  = wave,
    dominant_period       = confirmed_period,
    rhythm_strength_qp    = NA,        # fill from chi-sq output
    mesor                 = cosinor_params$MESOR,
    amplitude             = cosinor_params$amplitude,
    acrophase_hour        = cosinor_params$acrophase,
    mean_n_changepoints   = mean(sapply(cp_results, function(d)
                                  length(d$change_points)), na.rm = TRUE),
    sd_cp_timing          = sd(all_cps, na.rm = TRUE),  # regularity
    arrhythmic            = is.na(confirmed_period)
  )
}
```

---

## Stage 5: Longitudinal / Growth Curve Modeling

With features extracted at ages 12, 14, and 16 for all participants, fit
multilevel growth curve models to examine developmental trajectories.

```r
library(lme4)
library(lmerTest)

# Example: Does acrophase (clock time of HR peak) shift later with age?
model_acrophase <- lmer(
  acrophase_hour ~ age + age_sq + sex + (1 + age | participant_id),
  data = features_long %>% mutate(age_sq = age^2)
)
summary(model_acrophase)

# Example: Does rhythm strength decline across adolescence?
model_strength <- lmer(
  rhythm_strength_qp ~ age * sex + (1 + age | participant_id),
  data = features_long
)
summary(model_strength)
```

### Key Longitudinal Research Questions

| Question | Model |
|---|---|
| Does circadian phase delay from age 12→16? | `acrophase ~ age + (age | id)` |
| Do sex differences in rhythms emerge during puberty? | `feature ~ age * sex + (age | id)` |
| Are early arrhythmic participants still arrhythmic at 16? | `arrhythmic ~ wave + (1 | id)` (GLMM) |
| Does rhythm amplitude decrease with age? | `amplitude ~ age + (age | id)` |
| Does intra-day regularity (SD of CP timing) change? | `sd_cp_timing ~ age + (age | id)` |

---

## Packages Required

Install all dependencies before beginning:

```r
# CRAN packages
install.packages(c(
  "tidyverse",
  "lubridate",
  "imputeTS",
  "depmixS4",
  "cosinor",
  "cosinor2",
  "Rbeast",
  "lme4",
  "lmerTest",
  "changepoint",
  "bcp"
))

# GitHub packages
devtools::install_github("hadjamar/BayesSpec")  # AutoNOM

# Optional: Python CyHMM via reticulate
# pip install numpy scipy (in your Python env)
```

---

## File / Directory Structure

```
project/
├── data/
│   ├── raw/
│   │   ├── wave1/          # one CSV per participant, age 12
│   │   ├── wave2/          # age 14
│   │   └── wave3/          # age 16
│   └── processed/
│       ├── imputed/        # output of Stage 1
│       └── features/       # output of Stage 4
├── R/
│   ├── 01_preprocess.R
│   ├── 02_periodicity.R
│   ├── 03_changepoints.R
│   ├── 04_features.R
│   └── 05_longitudinal.R
├── outputs/
│   ├── figures/
│   └── tables/
└── handoff_biological_rhythms.md   # this document
```

---

## Known Issues and Decisions to Resolve

| Issue | Status | Recommendation |
|---|---|---|
| CyHMM R implementation is approximate | Open | Use Python via `reticulate` for fidelity |
| AutoNOM `BayesSpec` may not install on all systems | Open | Use `Rbeast` as fallback |
| Fitbit HR data noisier than E4 used in paper | Known | Accept attenuated rhythm strength estimates; note in limitations |
| Non-wear vs. true missing not always distinguishable | Open | Use Fitbit's own wear-detection output if available |
| kmax for AutoNOM selected empirically (k=4 in paper) | Open | Re-derive using MAPE elbow on your data; may differ for adolescents |
| Lomb-Scargle not in original paper but recommended for missing data | Enhancement | Add as 4th periodicity method given your missingness pattern |

---

## References

- Yan, R., & Doryab, A. (2022). Towards a Computational Framework for Automated
  Discovery and Modeling of Biological Rhythms from Wearable Data Streams.
  *IntelliSys 2021*, LNNS 296, 643–661. https://doi.org/10.1007/978-3-030-82199-9_44
- Hadj-Amar, B., et al. (2019). Bayesian Model Search for Nonstationary Periodic
  Time Series. *JASA*. [AutoNOM source]
- Pierson, E., Althoff, T., & Leskovec, J. (2018). Modeling Individual Cyclic Variation
  in Human Behavior. *WWW 2018*. [CyHMM source]
- Sokolove, P., & Bushell, W. (1978). The chi square periodogram. *J. Theor. Biol.*
- Moritz, S., & Bartz-Beielstein, T. (2017). imputeTS. *R Journal*.
