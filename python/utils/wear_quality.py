"""
Non-wear detection utilities.

TODO: sleep_stages — full implementation requires raw per-minute sleep stage
data which is not yet available. This stub flags nothing as non-wear.
Search for '# TODO: sleep_stages' to find all related stubs.
"""
import pandas as pd


def flag_nonwear(df: pd.DataFrame) -> pd.Series:
    """
    Return a boolean Series (True = non-wear slot) aligned with df.index.

    TODO: sleep_stages — replace with sustained flat-HR + zero-steps detection
          once raw per-minute data is available.
    """
    return pd.Series(False, index=df.index)
