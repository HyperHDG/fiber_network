# fiber_network

Python tools for dense fiber networks as used by [HyperHDG](https://github.com/HyperHDG/HyperHDG):
generate or convert a network into a `.geo.h5` domain, look at it, and compute its statistics.

## How to install

Python 3.10 or newer. The packages (numpy, scipy, pandas, h5py, pyyaml, matplotlib, and
pytest for the tests) are listed in `requirements.txt`:

```
python -m venv .venv && source .venv/bin/activate     # optional
pip install -r requirements.txt
```

`net_vis.py` and `net_vis_coarse.py` render with ParaView and run under its batch interpreter
`pvbatch` (or `pvpython` for interactive windows). Both come with every ParaView install:
`pacman -S paraview`, `apt install paraview`, `brew install --cask paraview`, or the binaries
from paraview.org.

The scripts run from a checkout, `python fiber_network/make_geo.py ...`, nothing needs to be
installed. To import the package from elsewhere, for an extension or from HyperHDG, install it
in place:

```
pip install -e .          # then: from fiber_network.make_geo import Network, make_parser, main
```

which also puts `make_geo` and `net_stats` on the PATH.

## Simplest example

```
python fiber_network/make_geo.py --mikado 200 -o mikado.geo.h5    # random fiber network on the unit square
python fiber_network/net_stats.py mikado.geo.h5                   # its statistics as YAML
pvbatch fiber_network/net_vis.py mikado.geo.h5 --warp-by none --arrows 0 -o mikado.png
```

## The domain file

A domain is one HDF5 file, `name.geo.h5`, with the group `/domain`:

| dataset | shape | content |
|---|---|---|
| `points` | (n_nodes, 3) | node coordinates |
| `edges` | (n_edges, 2) | node indices of each beam |
| `properties` | (n_edges, 17) | per-edge material data: mass, six stiffnesses, two normals, two widths, fiber id and edge id within the fiber, ... |
| `types_points`, `types_faces` | (n_nodes,), (n_edges, 2) | type of each node and of each edge end, 0 inside. Ternary by default: 1 or 2 on the x faces, 3 or 6 on y, 9 or 18 on z, summed on corners, as HyperHDG's cubic node descriptor; or a bit mask over the six beam dofs with `--dirichlet-bitmask` |

Attributes on `/domain` carry the sheet size and units. A `/VTKHDF` group holds a view of the
same arrays, so ParaView opens the file directly. `make_geo.py` can also write the legacy
ASCII `.geo` format that HyperHDG read before.

## Scripts

- `make_geo.py`: make a domain. Generators: `--grid N`, `--hex N`, `--mikado MASS [R]`;
  or convert a Morgan csv network (`-i folder` with `nodes.csv`, `fibers.csv`,
  `connections.csv`). Processing, in this order: clamp to a box (`--clamp-xy`), merge
  coincident nodes, subdivide long edges (`--subdivide P`, no edge longer than 1/P of the
  sheet), property cutoffs, rescale, node types for the box faces (`--dirichlet xmin ymax`,
  HyperHDG's ternary convention; or `--dirichlet-bitmask xmin=0b111111 center=0b000111` with
  one bit per beam dof, the timowave convention). `-o x.geo` selects the legacy format.
- `net_stats.py`: descriptive statistics of one or more domains as YAML, one top-level key per
  file: sheet, fibers, bonds, z distribution, cuts, optional METIS partition sizes. Written for
  real paper networks; on synthetic grids the fiber and bond sections are mostly empty.
- `net_vis.py` (pvbatch): render and animate a domain or a VTKHDF result on it: color by a
  field, warp by displacement, arrows, tubes, views, single frames or videos.
- `net_vis_coarse.py` (pvbatch): coarse 2D histogram of a VTKHDF field over the sheet, one
  image or an animation over time.
- `degree_distribution.py`: node degree distribution, optionally resolved in space.
- `distribution_length_stiffness.py`: 2D histogram of edge length against one property column.
- `dirichlet_spacing_vis.py`: the Dirichlet-node spacing outliers along the clamped faces.

Every script prints its options with `--help`.

## Extending make_geo

A generator that does not belong in this repo is added without touching it. `make_geo.py`
exposes three things: `Network`, whose `generate_from_args(args, parser)` method turns parsed
options into a network; `make_parser()`, the command line; and `main(argv, network_cls, parser)`,
the pipeline. A subclass adds its options to the parser, overrides `generate_from_args` to
handle them and defers to `super()` otherwise, and passes itself and the parser to `main`.
Processing, node types and output are inherited. `examples/custom_generator.py` does this for
a ladder graph:

```
python examples/custom_generator.py --ladder 8 --subdivide 4 -o ladder.geo.h5
```

The hook returns `True` for a generated network and `False` for one read from files, which is
what decides whether floating and small components are pruned.

## Examples

`examples/` has three scripts that write into `examples/output/`:

```
examples/generate.sh     # grid, honeycomb, mikado networks; subdivision and clamping; legacy .geo
examples/analyze.sh      # net_stats, degree distribution, histograms, Dirichlet spacing
examples/visualize.sh    # pvbatch renders
examples/custom_generator.py --ladder 8 -o ladder.geo.h5    # a generator added from outside
```

## Tests

```
python -m pytest
```

The tests generate small networks and run every script on them. The `net_vis.py` test is
skipped without `pvbatch`. The same suite runs on GitHub Actions for every push and pull
request (`.github/workflows/tests.yml`), without ParaView.

## Preconditioner

`gortz_hellman_malqvist_22.py`, the two-level domain decomposition preconditioner, still lives
here and is imported by HyperHDG as `fiber_network.precond`. It is moving to its own
repository.
