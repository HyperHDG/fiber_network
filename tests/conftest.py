import os
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "fiber_network"
sys.path.insert(0, str(SCRIPTS))   # the scripts import each other as top-level modules
sys.path.insert(0, str(ROOT / "examples"))

os.environ.setdefault("MPLBACKEND", "Agg")


def run(script, *args, **kw):
    """Run fiber_network/<script> with the test interpreter (or exe=pvbatch); returns CompletedProcess."""
    exe = kw.pop("exe", sys.executable)
    return subprocess.run([exe, str(SCRIPTS / script), *map(str, args)],
                          capture_output=True, text=True, check=True, **kw)


@pytest.fixture(scope="session")
def grid_h5(tmp_path_factory):
    out = tmp_path_factory.mktemp("grid") / "grid.geo.h5"
    run("make_geo.py", "--grid", 6, "-o", out)
    return out


@pytest.fixture(scope="session")
def mikado_h5(tmp_path_factory):
    out = tmp_path_factory.mktemp("mikado") / "mikado.geo.h5"
    run("make_geo.py", "--mikado", 40, "--seed", 1, "--subdivide", 4, "-o", out)
    return out
