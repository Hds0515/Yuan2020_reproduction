"""One-click v4 pipeline with strict output and protected-input isolation."""

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
    subprocess.run(
        [sys.executable, str(ROOT / script), *arguments],
        cwd=ROOT,
        check=True,
        timeout=timeout,
    )


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


def protected_paths() -> list[Path]:
    paths = sorted((ROOT / "source_images").glob("*"))
    paths += [
        ROOT / "identification" / "frozen_parameters_v2.json",
        ROOT / "identification_v3" / "frozen_parameters_v3.json",
    ]
    for pattern in ("*_raw_digitized.csv", "*_cleaned_1s.csv"):
        paths += sorted((ROOT / "digitization").glob(pattern))
    return [path for path in paths if path.is_file()]


def protected_manifest() -> dict[str, str]:
    return {path.relative_to(ROOT).as_posix(): digest(path) for path in protected_paths()}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, default=ROOT / "outputs_v4")
    parser.add_argument("--clean-output", action="store_true")
    parser.add_argument("--with-comsol", action="store_true")
    parser.add_argument("--comsol-timeout-s", type=int, default=3600)
    parser.add_argument("--safe-temperature-C", type=float)
    parser.add_argument(
        "--finalize-existing",
        action="store_true",
        help="Validate and manifest an already generated v4 bundle without rerunning experiments.",
    )
    args = parser.parse_args()
    output = args.output_dir.resolve()
    protected_before = protected_manifest()
    if args.finalize_existing:
        if not output.is_dir() or not any(output.iterdir()):
            raise FileNotFoundError(f"Existing output bundle not found: {output}")
    else:
        prepare_output(output, args.clean_output)
        run("digitization/digitize_figures.py", "--output-root", str(output))
        run("digitization/validate_digitization.py", "--output-root", str(output))
        run(
            "identification_v3/run_identification.py",
            "--root",
            str(ROOT),
            "--output-dir",
            str(output),
            "--digitized-dir",
            str(output / "digitization" / "regenerated"),
        )
        run(
            "experiments/run_observer_experiments.py",
            "--root",
            str(ROOT),
            "--output-dir",
            str(output),
        )
        control_arguments = ["--root", str(ROOT), "--output-dir", str(output)]
        if args.safe_temperature_C is not None:
            control_arguments += ["--safe-temperature-C", str(args.safe_temperature_C)]
        run("experiments/run_control_experiments.py", *control_arguments)
        subprocess.run([sys.executable, "-m", "pytest", "-q"], cwd=ROOT, check=True)

        if args.with_comsol:
            try:
                subprocess.run(
                    [str(ROOT / "comsol" / "run_comsol_windows.bat")],
                    cwd=ROOT / "comsol",
                    check=True,
                    timeout=args.comsol_timeout_s,
                    shell=True,
                )
            except (subprocess.CalledProcessError, subprocess.TimeoutExpired) as exception:
                print(f"COMSOL Stage 1 did not complete: {exception}", file=sys.stderr)

    run("finalize_results.py", "--output-dir", str(output))
    protected_after = protected_manifest()
    if protected_after != protected_before:
        changed = sorted(set(protected_before) | set(protected_after))
        changed = [name for name in changed if protected_before.get(name) != protected_after.get(name)]
        raise RuntimeError(f"Protected inputs changed during v4 run: {changed}")

    input_paths = [ROOT / "source_images" / f"Fig{number}.jpg" for number in (14, 15, 16)] + [
        ROOT / "digitization" / "axis_calibration.json",
        ROOT / "configs" / "reproduction_config.yaml",
        ROOT / "identification_v3" / "frozen_parameters_v3.json",
    ]
    try:
        git_commit = subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True
        ).strip()
    except (subprocess.CalledProcessError, FileNotFoundError):
        git_commit = None
    packages = {
        name: version(name)
        for name in ("numpy", "pandas", "scipy", "matplotlib", "Pillow", "PyYAML", "pytest")
    }
    provenance = {
        "python": sys.version,
        "operating_system": platform.platform(),
        "git_commit_at_run_start": git_commit,
        "packages": packages,
        "config_sha256": digest(ROOT / "configs" / "reproduction_config.yaml"),
        "protected_inputs_unchanged": True,
        "inputs": {
            str(path.relative_to(ROOT)).replace(os.sep, "/"): digest(path) for path in input_paths
        },
    }
    (output / "run_provenance.json").write_text(
        json.dumps(provenance, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    output_hashes = {
        path.relative_to(output).as_posix(): digest(path)
        for path in sorted(output.rglob("*"))
        if path.is_file() and path.name != "output_manifest_sha256.json"
    }
    (output / "output_manifest_sha256.json").write_text(
        json.dumps(output_hashes, indent=2) + "\n", encoding="utf-8"
    )
    print(f"Fresh v4 result bundle: {output}")


if __name__ == "__main__":
    main()
