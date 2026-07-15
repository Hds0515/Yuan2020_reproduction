# V2 results (read-only reference)

The immutable V2 source and result state is preserved by Git tag `baseline_v2`.
This directory is intentionally a reference rather than a copied result bundle, so a v4 run cannot accidentally include V2 artifacts. The legacy tracked `outputs/` and `outputs_v3/` directories are not read or copied by `run_all.py`; all fresh results are generated directly under `outputs_v4/`.
