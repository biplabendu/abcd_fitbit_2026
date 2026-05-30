"""
Run the full Python preprocessing and clustering pipeline in order.
Each script is run as a subprocess so log files are written correctly.

Usage:
    python python/run_pipeline.py
"""
import subprocess
import sys
from pathlib import Path

SCRIPTS = [
    "01_artifact_removal.py",
    "02_imputation.py",
    "03_valid_day_classification.py",
    "04_feature_extraction.py",
    "05_dimensionality_reduction.py",
    "06_clustering.py",
]

python     = sys.executable
script_dir = Path(__file__).parent

for script in SCRIPTS:
    print(f"\n{'=' * 60}\n  {script}\n{'=' * 60}")
    subprocess.run(
        [python, str(script_dir / script)],
        check=True,
    )

print("\nPipeline complete.")
