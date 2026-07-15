"""Five-node model, two-sensor EKF, and hotspot-constrained MPC.

The executable result files are in ../outputs.  This module is intentionally
kept separate from the frozen three-node baseline so later research extensions
do not change the reproduction baseline.
"""
from pathlib import Path
import json

ROOT = Path(__file__).resolve().parents[1]
print(json.loads((ROOT / "outputs" / "five_node_observer_mpc_metrics.json").read_text(encoding="utf-8")))
