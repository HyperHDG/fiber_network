#!/bin/bash
# Statistics and histograms of a domain. Run examples/generate.sh first.
set -eo pipefail
cd "$(dirname "$0")/.."
O=examples/output; D=$O/mikado200.geo.h5
export MPLBACKEND=${MPLBACKEND:-Agg}      # drop this line to get interactive windows

python fiber_network/net_stats.py $D $O/grid20.geo.h5 > $O/stats.yaml   # one top-level key per file
head -30 $O/stats.yaml

python fiber_network/degree_distribution.py $D                            # table on stdout
python fiber_network/degree_distribution.py $D --spatial --gridsize 8 --save $O/degree-spatial.png

# 2D histogram of edge length against property column 1 (first stiffness); column 0 is the mass
python fiber_network/distribution_length_stiffness.py -i $D -c 1 -b 20 -o $O/length-stiffness.png --no-show

# where the Dirichlet nodes on the clamped faces are spaced widest
python fiber_network/dirichlet_spacing_vis.py $D -o $O/dirichlet-spacing.png
ls -la $O/*.png $O/stats.yaml
