from pathlib import Path

import pandas as pd

from digitization.digitize_figures import digitize_figure


def test_fig15_is_regenerated_from_image(tmp_path: Path) -> None:
    metrics = digitize_figure("Fig15", tmp_path)
    output = tmp_path / "digitization" / "regenerated" / "Fig15_digitized.csv"
    data = pd.read_csv(output)
    assert list(data.columns) == ["time_s", "experiment_current_A", "paper_simulation_current_A"]
    assert len(data) == 1201
    assert metrics["curves"]["paper_simulation_current_A"]["valid_column_ratio"] > 0.95
    assert (tmp_path / "digitization" / "Fig15_pixel_overlay.png").exists()
