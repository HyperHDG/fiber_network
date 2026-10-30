"""make_geo.py: generators and processing steps through the class, the whole
pipeline through the command line."""
import re
import subprocess
import sys

import h5py
import numpy as np

import make_geo
from conftest import ROOT, run


def check_graph(net):
    e = net.edges
    assert e.min() >= 0 and e.max() < len(net.nodes)
    assert (e[:, 0] != e[:, 1]).all(), "self loop"
    assert len(np.unique(np.sort(e, axis=1), axis=0)) == len(e), "duplicate edge"


def edge_lengths(net):
    p = net.nodes[net.edges]
    return np.linalg.norm(p[:, 1] - p[:, 0], axis=1)


def test_grid_counts_and_extent():
    net = make_geo.Network()
    net.generate_grid(5, 4)
    assert net.nodes.shape == (20, 3)
    assert net.edges.shape == ((5 - 1) * 4 + 5 * (4 - 1), 2)
    assert net.nodes[:, :2].min() == 0 and net.nodes[:, :2].max() == 1
    check_graph(net)


def test_honeycomb_is_a_graph():
    net = make_geo.Network()
    net.generate_honeycomb(3, 3)
    assert len(net.nodes) > 0
    check_graph(net)


def test_mikado_seed_reproducible_and_inside_unit_square():
    a, b, c = make_geo.Network(), make_geo.Network(), make_geo.Network()
    a.generate_mikado(30, seed=2)
    b.generate_mikado(30, seed=2)
    c.generate_mikado(30, seed=3)
    np.testing.assert_array_equal(a.nodes, b.nodes)
    np.testing.assert_array_equal(a.edges, b.edges)
    assert not (a.nodes.shape == c.nodes.shape and np.array_equal(a.nodes, c.nodes))
    assert a.nodes[:, :2].min() >= 0 and a.nodes[:, :2].max() <= 1
    check_graph(a)


def test_synthetic_properties_mass_is_length():
    net = make_geo.Network()
    net.generate_grid(4, 4)
    net.generate_synthetic_properties(width=0.01)
    assert net.edgeProps.shape == (len(net.edges), 17)
    np.testing.assert_allclose(net.edgeProps[:, 0], edge_lengths(net))


def test_subdivide_bounds_edge_length_and_conserves_mass():
    net = make_geo.Network()
    net.generate_grid(3, 3)                      # edge length 1/2
    net.generate_synthetic_properties(width=0.01)
    mass = net.edgeProps[:, 0].sum()
    net.subdivide_long_edges(4)                  # no edge longer than 1/4
    assert edge_lengths(net).max() <= 0.25 * (1 + 1e-12)
    assert len(net.edges) == 2 * 12
    np.testing.assert_allclose(net.edgeProps[:, 0].sum(), mass)
    check_graph(net)


def test_rescale_bbox_and_scale():
    net = make_geo.Network()
    net.generate_mikado(30, seed=0)
    net.scale(3.0)
    assert net.nodes[:, :2].max() > 1
    net.rescale_bbox()
    ext = net.nodes[:, :2].max(0) - net.nodes[:, :2].min(0)
    np.testing.assert_allclose(ext.max(), 1.0)
    np.testing.assert_allclose(net.info["size"][:2], ext)


# command line: the fixtures in conftest.py ran make_geo.py --grid 6 and --mikado 40 --subdivide 4

def test_cli_geo_h5_layout_and_dirichlet_types(grid_h5):
    with h5py.File(grid_h5) as f:
        d = f["domain"]
        assert d["points"].shape == (36, 3)
        assert d["edges"].shape == (60, 2)
        assert d["properties"].shape == (60, 17)
        types_points = d["types_points"][:]
        types_faces = d["types_faces"][:]
        assert "VTKHDF" in f, "ParaView view missing"
    # default --dirichlet xmin xmax: ternary types 1 and 2, two columns of 6 nodes
    assert (types_points != 0).sum() == 12
    assert (types_points == 1).sum() == 6 and (types_points == 2).sum() == 6
    assert types_faces.shape == (60, 2)


def test_cli_dirichlet_ternary_corner(tmp_path):
    out = tmp_path / "g.geo.h5"
    run("make_geo.py", "--grid", 6, "--dirichlet", "xmin", "ymax", "-o", out)
    with h5py.File(out) as f:
        t = f["domain/types_points"][:]
    assert (t == 1).sum() == 5 and (t == 6).sum() == 5 and (t == 7).sum() == 1   # 1 + 2*3 on the corner
    assert (t != 0).sum() == 11


def test_cli_dirichlet_bitmask(tmp_path):
    out = tmp_path / "g.geo.h5"
    run("make_geo.py", "--grid", 6, "--dirichlet-bitmask", "xmin=0b111111", "ymax=0b000111", "-o", out)
    with h5py.File(out) as f:
        t = f["domain/types_points"][:]
    assert (t == 0b111111).sum() == 6 and (t == 0b000111).sum() == 5   # masks OR on the corner


def test_cli_legacy_geo(tmp_path):
    out = tmp_path / "grid.geo"
    run("make_geo.py", "--grid", 4, "-o", out)
    text = out.read_text()
    assert re.search(r"N_Points\s*=\s*16;", text)
    assert re.search(r"N_HyperEdges\s*=\s*24;", text)
    points = text.split("POINTS:")[1].split("HYPERNODES_OF_HYPEREDGES:")[0].split()
    assert len(points) == 16 * 3
    assert "HYPEREDGE_PROPERTIES: 17" in text


def test_cli_subdivide(mikado_h5):
    with h5py.File(mikado_h5) as f:
        points, edges = f["domain/points"][:], f["domain/edges"][:]
    ext = (points[:, :2].max(0) - points[:, :2].min(0)).max()
    lengths = np.linalg.norm(points[edges[:, 1]] - points[edges[:, 0]], axis=1)
    assert lengths.max() <= ext / 4 * (1 + 1e-12)


# node types through the class (no command line involved)

def test_compute_types_ternary_and_bitmask():
    net = make_geo.Network()
    net.generate_grid(4, 4)
    net.compute_types(0.02, faces=("ymin",))
    assert (net.types_points == 3).sum() == 4 and (net.types_points != 0).sum() == 4
    net.compute_types(0.02, faces=("xmin", "zmax"))          # flat network: no z faces
    assert set(np.unique(net.types_points)) == {0, 1}
    net.compute_types(0.02, bitmask=["xmin=0b11", "center=0b100"], center_radius=0.3)
    assert (net.types_points & 0b11 == 0b11).sum() == 4
    assert (net.types_points & 0b100).sum() > 0
    assert net.types_faces.shape == (len(net.edges), 2)


# extending the generator the way examples/custom_generator.py does

def test_extension_in_process(tmp_path):
    import custom_generator as cg
    out = tmp_path / "ladder.geo.h5"
    net = make_geo.main(["--ladder", "5", "-o", str(out)], network_cls=cg.LadderNetwork, parser=cg.parser)
    assert isinstance(net, cg.LadderNetwork)
    assert len(net.nodes) == 10 and len(net.edges) == 2 * 4 + 5
    with h5py.File(out) as f:
        assert f["domain/points"].shape == (10, 3)
        assert (f["domain/types_points"][:] != 0).sum() == 4    # default xmin xmax: both rail ends
    # the inherited generators and options keep working through the extended parser
    net = make_geo.main(["--grid", "3", "--subdivide", "4", "-o", str(out)], network_cls=cg.LadderNetwork, parser=cg.parser)
    assert len(net.edges) == 24


def test_extension_cli(tmp_path):
    out = tmp_path / "ladder.geo.h5"
    subprocess.run([sys.executable, str(ROOT / "examples" / "custom_generator.py"),
                    "--ladder", "6", "--dirichlet", "xmin", "-o", str(out)],
                   check=True, capture_output=True, text=True)
    with h5py.File(out) as f:
        t = f["domain/types_points"][:]
    assert len(t) == 12 and (t == 1).sum() == 2
