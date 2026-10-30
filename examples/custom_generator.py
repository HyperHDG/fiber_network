#!/usr/bin/env python3
"""A generator of your own, outside this repo.

Subclass Network, add your options to the parser from make_parser(), override
generate_from_args to handle them, and hand both to main(). Everything else, the
processing steps, node types and output formats, is inherited:

    python examples/custom_generator.py --ladder 8 -o ladder.geo.h5
    python examples/custom_generator.py --ladder 8 --subdivide 4 --dirichlet xmin -o ladder.geo.h5
    python examples/custom_generator.py --grid 5 -o grid.geo.h5     # the inherited generators still work

From a checkout the import below finds make_geo.py next to this file's parent; an
installed package (pip install -e .) is imported as fiber_network.make_geo instead.
"""
import pathlib
import sys

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "fiber_network"))
from make_geo import Network, make_parser, main   # noqa: E402


class LadderNetwork(Network):

  def generate_ladder(self, n, width=0.2):
    """Two rails along x on the unit square, n rungs between them."""
    x = np.linspace(0.0, 1.0, n)
    nodes = np.zeros((2 * n, 3))
    nodes[:n, 0] = x
    nodes[n:, 0] = x
    nodes[n:, 1] = width
    rails = [(i, i + 1) for i in range(n - 1)] + [(n + i, n + i + 1) for i in range(n - 1)]
    rungs = [(i, n + i) for i in range(n)]
    self.nodes = nodes
    self.edges = np.array(rails + rungs, dtype=np.int64)
    self.info = {"size": np.array([1.0, width, 0.0])}
    self.edgeProps = None

  def generate_from_args(self, args, parser):
    if args.ladder is None:
      return super().generate_from_args(args, parser)
    self.generate_ladder(args.ladder, width=args.ladder_width)
    self.generate_synthetic_properties(width=0.01)
    return True           # generated, not read: no component pruning


parser = make_parser()
parser.add_argument("--ladder", type=int, metavar="N", help="generate a ladder with N rungs")
parser.add_argument("--ladder-width", type=float, default=0.2, help="distance of the rails")

if __name__ == "__main__":
  main(network_cls=LadderNetwork, parser=parser)
