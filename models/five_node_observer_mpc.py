"""Compatibility entry point for the executable five-node research extension.

The implementation is intentionally split into auditable model, observer,
controller, and experiment modules.  Running this file executes both experiment
suites; it never reads a pre-existing metrics JSON as a substitute for work.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=ROOT)
    args = parser.parse_args(); root = args.root.resolve()
    for script in ("experiments/run_observer_experiments.py", "experiments/run_control_experiments.py"):
        subprocess.run([sys.executable, str(root / script), "--root", str(root)], cwd=root, check=True)


if __name__ == "__main__":
    main()
