#!/bin/bash
# Pictures with ParaView's pvbatch (offscreen; pvbatch for --show windows). Run examples/generate.sh first.
set -eo pipefail
cd "$(dirname "$0")/.."
O=examples/output; D=$O/mikado200-sub8.geo.h5

# the domain itself: no field to warp by, no arrows. Edges are drawn as hollow rectangular
# beams by default (--beams 1), with normals and widths from the property columns
# (--beams-cols, default 7,8,9:10,11,12:13:14); --beams-scale widens them for a thin network
pvbatch fiber_network/net_vis.py $D --warp-by none --arrows 0 -o $O/mikado.png
pvbatch fiber_network/net_vis.py $D --warp-by none --arrows 0 --view iso --beams-scale 3 -o $O/mikado-iso.png

# tubes of a fixed radius instead of beams, and a filled surface for networks finer than pixels
pvbatch fiber_network/net_vis.py $D --warp-by none --arrows 0 --beams 0 -r 0.003 --tubes-sides 8 -o $O/mikado-tubes.png
pvbatch fiber_network/net_vis.py $D --warp-by none --arrows 0 --surface -o $O/mikado-surface.png

# a solver result on the domain (a VTKHDF file with point fields, e.g. from HyperHDG) colors
# and warps by a field; "values:0" is component 0 of the array "values":
#   pvbatch fiber_network/net_vis.py result.vtkhdf --color-by values:2 --warp-by values -r 20 -o result.png
#   pvbatch fiber_network/net_vis.py wave.vtkhdf --frames 0,0.3,0.6 -o wave.mp4
# and a coarse 2D histogram of that field over the sheet:
#   pvbatch fiber_network/net_vis_coarse.py wave.vtkhdf --color-by values:8 --bins 60 -o hist.png
ls -la $O/*.png
