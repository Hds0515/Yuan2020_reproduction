"""One-click reproducibility entry point starting from the protected JPG files."""

from __future__ import annotations

import argparse
import hashlib
from importlib.metadata import version
import json
import os
from pathlib import Path
import platform
import shutil
import subprocess
import sys

ROOT = Path(__file__).resolve().parent


def run(script: str, *arguments: str, timeout: int | None = None) -> None:
    command = [sys.executable, str(ROOT / script), *arguments]
    subprocess.run(command, cwd=ROOT, check=True, timeout=timeout)


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def prepare_output(path: Path, clean: bool) -> None:
    resolved = path.resolve()
    if ROOT not in resolved.parents or resolved == ROOT:
        raise ValueError("Output directory must be a child of the project root")
    if resolved.exists() and any(resolved.iterdir()):
        if not clean:
            raise FileExistsError(f"Output directory is not empty: {resolved}; pass --clean-output")
        shutil.rmtree(resolved)
    resolved.mkdir(parents=True, exist_ok=True)


def clear_generated_workspace() -> None:
    """Remove only v3-generated artifacts; never touch protected V2 inputs."""
    generated_directories = [
        ROOT / "digitization" / "regenerated",
        ROOT / "outputs" / "digitization",
        ROOT / "outputs" / "three_node",
        ROOT / "outputs" / "five_node",
        ROOT / "outputs" / "observer",
        ROOT / "outputs" / "control",
    ]
    for path in generated_directories:
        if path.exists():
            shutil.rmtree(path)
    for path in (ROOT / "outputs" / "final_summary.json", ROOT / "outputs" / "final_summary.csv"):
        if path.exists():
            path.unlink()
    identification = ROOT / "identification_v3"
    for path in identification.iterdir():
        if path.is_file() and path.name != "run_identification.py":
            path.unlink()


def collect(output: Path) -> None:
    for source, destination in (
        (ROOT / "digitization" / "regenerated", output / "digitization" / "regenerated"),
        (ROOT / "outputs", output / "outputs"),
        (ROOT / "docs", output / "docs"),
    ):
        shutil.copytree(source, destination, dirs_exist_ok=True)
    identification_destination = output / "identification_v3"; identification_destination.mkdir(parents=True, exist_ok=True)
    for path in (ROOT / "identification_v3").iterdir():
        if path.is_file() and path.suffix != ".py":
            shutil.copy2(path, identification_destination / path.name)
    shutil.copy2(ROOT / "configs" / "reproduction_config.yaml", output / "reproduction_config.yaml")
    shutil.copy2(ROOT / "audit" / "final_manifest_sha256.json", output / "final_manifest_sha256.json")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, default=ROOT / "outputs_v3")
    parser.add_argument("--clean-output", action="store_true")
    parser.add_argument("--with-comsol", action="store_true")
    parser.add_argument("--comsol-timeout-s", type=int, default=1800)
    args = parser.parse_args(); output = args.output_dir.resolve()
    prepare_output(output, args.clean_output)
    clear_generated_workspace()

    run("digitization/digitize_figures.py", "--output-root", str(ROOT))
    # V2 comparison is diagnostic only and runs strictly after image extraction.
    run("digitization/validate_digitization.py", "--output-root", str(ROOT))
    run("identification_v3/run_identification.py", "--root", str(ROOT))
    run("experiments/run_observer_experiments.py", "--root", str(ROOT))
    run("experiments/run_control_experiments.py", "--root", str(ROOT))
    subprocess.run([sys.executable, "-m", "pytest", "-q"], cwd=ROOT, check=True)
    if args.with_comsol:
        try:
            subprocess.run([str(ROOT / "comsol" / "run_comsol_windows.bat")], cwd=ROOT / "comsol",
                           check=True, timeout=args.comsol_timeout_s, shell=True)
        except (subprocess.CalledProcessError, subprocess.TimeoutExpired) as exception:
            print(f"COMSOL did not complete: {exception}", file=sys.stderr)
    run("finalize_results.py")
    run("audit/generate_final_manifest.py")
    collect(output)

    input_paths = [ROOT / "source_images" / f"Fig{number}.jpg" for number in (14, 15, 16)] + [
        ROOT / "digitization" / "axis_calibration.json", ROOT / "configs" / "reproduction_config.yaml"]
    try:
        git_commit = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
    except (subprocess.CalledProcessError, FileNotFoundError):
        git_commit = None
    packages = {name: version(name) for name in ("numpy", "pandas", "scipy", "matplotlib", "Pillow", "PyYAML", "pytest")}
    provenance = {
        "python": sys.version,
        "operating_system": platform.platform(),
        "git_commit_at_run_start": git_commit,
        "packages": packages,
        "config_sha256": digest(ROOT / "configs" / "reproduction_config.yaml"),
        "inputs": {str(path.relative_to(ROOT)).replace(os.sep, "/"): digest(path) for path in input_paths},
    }
    (output / "run_provenance.json").write_text(json.dumps(provenance, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    output_hashes = {}
    for path in sorted(output.rglob("*")):
        if path.is_file() and path.name not in {"run_provenance.json", "output_manifest_sha256.json"}:
            output_hashes[path.relative_to(output).as_posix()] = digest(path)
    (output / "output_manifest_sha256.json").write_text(json.dumps(output_hashes, indent=2) + "\n", encoding="utf-8")
    print(f"Fresh result bundle: {output}")


if __name__ == "__main__":
    main()
