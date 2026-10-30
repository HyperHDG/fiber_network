#!/bin/bash
# Generate domains with make_geo.py. Output goes to examples/output/.
set -eo pipefail
cd "$(dirname "$0")/.."
G="python fiber_network/make_geo.py"
mkdir -p examples/output; O=examples/output

$G --grid 20          -o $O/grid20.geo.h5        # 20 x 20 grid graph on the unit square
$G --hex 12 8         -o $O/hex12x8.geo.h5       # honeycomb, 12 x 8 cells
$G --mikado 200       -o $O/mikado200.geo.h5     # random fibers of total length 200, fiber length 0.05
$G --mikado 200 0.1 --seed 3 -o $O/mikado200-long.geo.h5   # longer fibers, another seed

# processing steps, same for every input:
#   --subdivide 8   no edge longer than 1/8 of the sheet (for 8 x 8 domain decomposition)
#   --rescale-bbox  xy bounding box to 1 x 1
#   --dirichlet     which box faces get a node type, HyperHDG's ternary convention
#                   (xmin 1, xmax 2, ymin 3, ymax 6, summed on corners); default xmin xmax
#   --dirichlet-bitmask FACE=MASK  instead: one bit per beam dof (timowave), also center=MASK
$G --mikado 200 --subdivide 8 --dirichlet ymin ymax -o $O/mikado200-sub8.geo.h5
$G --mikado 200 --dirichlet-bitmask xmin=0b111111 center=0b000111 -o $O/mikado200-bitmask.geo.h5

# legacy ASCII .geo (the format HyperHDG read before .geo.h5)
$G --grid 6 -o $O/grid6.geo

# a Morgan csv network (nodes.csv, fibers.csv, connections.csv in one folder) converts with
#   $G -i path/to/folder --clamp-xy 0.5 --rescale-bbox -o net.geo.h5
ls -la $O
