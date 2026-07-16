"""Generate a stable SHA-256 manifest, excluding volatile environments/artifacts."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "audit" / "final_manifest_sha256.json"
EXCLUDED_PARTS = {".git", ".venv", "__pycache__", ".pytest_cache", "outputs_v3"}


def main() -> None:
    records = []
    for path in sorted(ROOT.rglob("*")):
        relative = path.relative_to(ROOT)
        if not path.is_file() or any(part in EXCLUDED_PARTS for part in relative.parts):
            continue
        if path == OUTPUT or path.suffix in {".class", ".mph"} or path.name.endswith(".class.status"):
            continue
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        records.append({"path": relative.as_posix(), "bytes": path.stat().st_size, "sha256": digest})
    payload = {"algorithm": "SHA-256", "self_excluded": True, "file_count": len(records), "files": records}
    OUTPUT.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(f"manifest files: {len(records)}")


if __name__ == "__main__":
    main()
