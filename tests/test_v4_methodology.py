from pathlib import Path

import yaml

from experiments.run_control_experiments import _safe_metrics


ROOT = Path(__file__).resolve().parents[1]


def test_reference_and_safe_limit_are_not_aliased() -> None:
    config = yaml.safe_load((ROOT / "configs" / "reproduction_config.yaml").read_text(encoding="utf-8"))
    assert config["model"]["reference_temperature_C"] == 55.0
    assert config["safety"]["safe_temperature_C"] is None
    metrics = _safe_metrics(hotspot=[54.0, 56.0], safe_C=None, dt_s=1.0)  # type: ignore[arg-type]
    assert metrics["safety_conclusion_permitted"] is False
    assert metrics["time_above_safe_limit_s"] is None


def test_publication_bootstrap_count_is_at_least_200() -> None:
    config = yaml.safe_load((ROOT / "configs" / "reproduction_config.yaml").read_text(encoding="utf-8"))
    assert config["identification"]["bootstrap_count"] >= 200


def test_v4_runner_does_not_collect_legacy_outputs() -> None:
    source = (ROOT / "run_all.py").read_text(encoding="utf-8")
    assert "copytree" not in source
    assert 'ROOT / "outputs_v4"' in source
