"""
Circadian rhythm analysis utilities.

Public API
----------
fit_cosinor(times, values)          -> (mesor, amplitude, acrophase_hours)
interdaily_stability(values_2d)     -> float
intradaily_variability(values_2d)   -> float
l5_m10(mean_profile, slot_hours)    -> dict
"""
import numpy as np


def fit_cosinor(
    times: np.ndarray,
    values: np.ndarray,
    period: float = 24.0,
) -> tuple[float, float, float]:
    """
    Linearised cosinor regression:  y = M + β·cos(2πt/T) + γ·sin(2πt/T)

    Parameters
    ----------
    times   : 1-D array of time-of-day values (hours, same length as values).
    values  : 1-D array of signal values (NaN-safe).
    period  : rhythm period in hours (default 24).

    Returns
    -------
    (mesor, amplitude, acrophase_hours)
        acrophase is in [0, period).  Returns (nan, nan, nan) when < 3 valid
        observations or when the least-squares fit fails.
    """
    valid = ~np.isnan(values)
    if valid.sum() < 3:
        return np.nan, np.nan, np.nan

    t, y = times[valid], values[valid]
    omega = 2 * np.pi / period
    X = np.column_stack([np.ones_like(t), np.cos(omega * t), np.sin(omega * t)])

    try:
        coeffs, _, _, _ = np.linalg.lstsq(X, y, rcond=None)
    except np.linalg.LinAlgError:
        return np.nan, np.nan, np.nan

    M, beta, gamma = coeffs
    amplitude = float(np.sqrt(beta**2 + gamma**2))
    acrophase_rad = np.arctan2(gamma, beta)   # φ = atan2(γ,β) for model y=M+β·cos+γ·sin
    acrophase_hours = float(acrophase_rad * period / (2 * np.pi) % period)

    return float(M), amplitude, acrophase_hours


def interdaily_stability(values_2d: np.ndarray) -> float:
    """
    IS = (n · Σ_h(x̄_h − x̄)²) / (p · Σ_i(x_i − x̄)²)

    Parameters
    ----------
    values_2d : shape (n_days, n_slots).  NaN values are ignored.

    Returns float in [0, 1], or NaN when insufficient data.
    """
    n_days, p = values_2d.shape
    flat = values_2d.flatten()
    n = int(np.sum(~np.isnan(flat)))
    if n < p:
        return np.nan

    grand_mean   = float(np.nanmean(flat))
    hourly_means = np.nanmean(values_2d, axis=0)   # (p,)

    numerator   = n * float(np.nansum((hourly_means - grand_mean) ** 2))
    denominator = p * float(np.nansum((flat - grand_mean) ** 2))

    return float(numerator / denominator) if denominator > 0 else np.nan


def intradaily_variability(values_2d: np.ndarray) -> float:
    """
    IV = (n · Σ(x_{i+1} − x_i)²) / ((n−1) · Σ(x_i − x̄)²)

    Computed on the full flattened time series (consecutive days).
    NaN transitions are skipped in the difference term.

    Returns float (typically 0–2), or NaN when insufficient data.
    """
    series = values_2d.flatten()
    valid_mask = ~np.isnan(series)
    n = int(valid_mask.sum())
    if n < 4:
        return np.nan

    grand_mean = float(np.nanmean(series))
    diffs = np.diff(series)
    sq_diffs = np.where(np.isnan(diffs), np.nan, diffs ** 2)

    numerator   = n * float(np.nansum(sq_diffs))
    denominator = (n - 1) * float(np.nansum((series - grand_mean) ** 2))

    return float(numerator / denominator) if denominator > 0 else np.nan


def l5_m10(
    mean_profile: np.ndarray,
    slot_hours: list,
    l5_slots: int = 3,
    m10_slots: int = 5,
) -> dict:
    """
    Compute L5 and M10 from a mean diurnal activity profile.

    Parameters
    ----------
    mean_profile : 1-D array of mean steps per slot (length = n_slots).
    slot_hours   : list of start hours for each slot (same length).
    l5_slots     : rolling window for L5 (default 3 ≈ 6 h at 2-hr resolution).
    m10_slots    : rolling window for M10 (default 5 = 10 h).

    Returns
    -------
    dict with keys: l5, l5_onset, m10, m10_onset, relative_amplitude.
    All values are NaN when the profile has too many missing slots.
    """
    p = len(mean_profile)
    nan_result = dict(l5=np.nan, l5_onset=np.nan, m10=np.nan,
                      m10_onset=np.nan, relative_amplitude=np.nan)

    max_window = max(l5_slots, m10_slots)
    if np.isnan(mean_profile).sum() > p - max_window:
        return nan_result

    # Double the profile for circular (wrap-around) window search
    ext = np.concatenate([mean_profile, mean_profile])

    def _rolling_mean(arr, w):
        return np.array([np.nanmean(arr[i : i + w]) for i in range(len(arr) - w + 1)])

    l5_means  = _rolling_mean(ext, l5_slots)[:p]
    m10_means = _rolling_mean(ext, m10_slots)[:p]

    l5_idx  = int(np.nanargmin(l5_means))
    m10_idx = int(np.nanargmax(m10_means))

    l5_val  = float(l5_means[l5_idx])
    m10_val = float(m10_means[m10_idx])
    ra = (m10_val - l5_val) / (m10_val + l5_val) if (m10_val + l5_val) > 0 else np.nan

    return dict(
        l5=l5_val,
        l5_onset=float(slot_hours[l5_idx]),
        m10=m10_val,
        m10_onset=float(slot_hours[m10_idx]),
        relative_amplitude=ra,
    )
