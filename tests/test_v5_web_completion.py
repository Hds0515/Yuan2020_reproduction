from pathlib import Path
import hashlib,json,yaml
ROOT=Path(__file__).resolve().parents[1]
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def test_primary_hash_and_runtime_records():
 c=yaml.safe_load((ROOT/"configs/reproduction_config.yaml").read_text())
 p=ROOT/c["primary_model"]["path"];assert digest(p)==c["primary_model"]["sha256"]
 for f in (ROOT/"outputs_v5").glob("*/model_runtime_*.json"):
  assert json.loads(f.read_text())["primary_model_sha256"]==c["primary_model"]["sha256"]
def test_negative_results_are_preserved():
 s=json.loads((ROOT/"outputs_v5/final_summary.json").read_text())
 assert s["distributed_model"]["five_node_external_validation_passed"] is False
 assert s["observer"]["MM_EKF_core_innovation_supported"] is False
 assert s["control"]["equal_energy_hotspot_advantage_proven"] is False
 assert s["comsol"]["web_environment_executed_comsol"] is False
def test_55C_is_reference_not_safety():
 c=yaml.safe_load((ROOT/"configs/reproduction_config.yaml").read_text())
 assert c["safety"]["safe_temperature_C"] is None
 assert "55°C仅为参考温度" in (ROOT/"docs_v5/control_fair_comparison.md").read_text()
