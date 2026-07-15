"""Single-source primary-model registry and hash verification for v5."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

import yaml


def sha256_file(path: str | Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load_config(root: str | Path) -> dict[str, Any]:
    root = Path(root)
    return yaml.safe_load((root / "configs" / "reproduction_config.yaml").read_text(encoding="utf-8"))


def primary_model_path(root: str | Path, config: dict[str, Any] | None = None) -> Path:
    root = Path(root).resolve()
    cfg = config or load_config(root)
    value = cfg.get("primary_model", {}).get("path")
    if not value:
        raise KeyError("configs/reproduction_config.yaml is missing primary_model.path")
    path = (root / value).resolve()
    if root not in path.parents:
        raise ValueError("primary_model.path must be inside the project root")
    if not path.is_file():
        raise FileNotFoundError(path)
    return path


def load_primary_model(root: str | Path) -> tuple[dict[str, Any], Path, str]:
    root = Path(root).resolve()
    config = load_config(root)
    path = primary_model_path(root, config)
    actual = sha256_file(path)
    expected = config.get("primary_model", {}).get("sha256")
    if expected and expected != actual:
        raise RuntimeError(
            f"Primary model hash mismatch: config={expected}, actual={actual}, path={path}"
        )
    data = json.loads(path.read_text(encoding="utf-8"))
    declared = data.get("file_sha256")
    if declared and declared != actual:
        raise RuntimeError(f"Primary model self-declared hash mismatch: {declared} != {actual}")
    return data, path, actual


def write_model_runtime_record(root: str | Path, output_dir: str | Path, consumer: str) -> Path:
    data, path, digest = load_primary_model(root)
    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    record = {
        "consumer": consumer,
        "primary_model_path": path.relative_to(Path(root).resolve()).as_posix(),
        "primary_model_sha256": digest,
        "heat_model": data["heat_model"],
        "model_name": data.get("model"),
    }
    target = output / f"model_runtime_{consumer}.json"
    target.write_text(json.dumps(record, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return target
