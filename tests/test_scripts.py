"""The analysis and visualization scripts, run as a user would run them."""
import re
import shutil

import pytest
import yaml

from conftest import run

PVBATCH = shutil.which("pvbatch")


def test_net_stats_yaml(grid_h5, mikado_h5):
    out = run("net_stats.py", grid_h5, mikado_h5).stdout
    doc = yaml.safe_load(out)
    assert set(doc) == {"grid", "mikado"}, "one top-level key per input file, named after it"
    sheet = doc["grid"]["sheet"]
    assert sheet["nodes"] == 36 and sheet["edges"] == 60
    assert {"sheet", "fibers", "bonds", "z"} <= set(doc["grid"])


def test_degree_distribution(grid_h5):
    out = run("degree_distribution.py", grid_h5).stdout
    assert re.search(r"nodes:\s+36", out)
    assert re.search(r"edges:\s+60", out)
    assert "max=4" in out           # interior grid nodes


def test_distribution_length_stiffness(grid_h5, tmp_path):
    png = tmp_path / "hist.png"
    run("distribution_length_stiffness.py", "-i", grid_h5, "-c", 1, "-o", png, "--no-show")
    assert png.stat().st_size > 0


def test_dirichlet_spacing_vis(mikado_h5, tmp_path):
    png = tmp_path / "spacing.png"
    run("dirichlet_spacing_vis.py", mikado_h5, "-o", png)
    assert png.stat().st_size > 0


@pytest.mark.skipif(PVBATCH is None, reason="needs ParaView's pvbatch")
def test_net_vis_renders(mikado_h5, tmp_path):
    png = tmp_path / "net.png"
    run("net_vis.py", mikado_h5, "--warp-by", "none", "--arrows", 0, "-o", png, exe=PVBATCH)
    assert png.stat().st_size > 0
