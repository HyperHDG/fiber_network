#!/usr/bin/env python

"""
net_stats.py -- descriptive statistics of a fibre network, as YAML.

Reads one or more .geo.h5 domains (make_geo2.py format: /domain/points, /edges,
/properties with the 17 timowave columns) and prints, per file, a YAML document
with the quantities needed to compare a generated network against the morgan data
and against Goertz et al. 2025 (Comput. Mech. 75:1-13, Sec. 3):

  sheet     bounding box, thickness, grammage and density (from the mass column)
  fibers    count, arc length, shape factor (= projected / arc length, Fig. 2),
            in-plane orientation and the fitted eta of the cos-1 PDF (eq. 4),
            out-of-plane angle, z-extent, segment length, cross-section columns
  bonds     count, length, |dz|, inclination, bonds per fibre pair, spacing between
            bonds of the same pair, bonds per node, bond columns
  z         node / fibre-length profile over the thickness, and edges cut by
            mid planes normal to x, y, z per unit plane area (a crude z-connectivity
            measure; the anisotropy ratio cut_z / cut_x is the headline number)
  partition (--metis K) METIS k-way partitions of the (unweighted) graph: edge cut and,
            per axis, the fraction of the node-coordinate variance explained by the
            parts (r2). r2_z near 1 means the parts are stacked in z, i.e. the network
            is easiest to cut through the thickness: too little z-connectivity.
            gpmetis is taken from $GPMETIS or PATH.

Units are those of the file (morgan / make_geo2: kg, um, s). Distributions are
reported as percentiles p5/p50/p95 (plus mean where it is meaningful).

Usage:
    net_stats.py golden/net1-sca.geo.h5
    net_stats.py a.geo.h5 b.geo.h5 > stats.yaml
    net_stats.py golden/net1-sca.geo.h5 --metis 2 3 4

Author
------
Joseph Holten, KIT, 2026.
"""

import argparse
import os
import shutil
import subprocess
import sys
import tempfile
import time

import h5py
import numpy as np
import scipy.sparse as sp
import yaml


def tprint(*args, **kwargs):
  print(f"[{time.strftime('%H:%M:%S')}]", *args, file=sys.stderr, **kwargs)


def num(x, sig=4):
  """numpy scalar -> plain python number rounded to `sig` significant digits."""
  x = float(x)
  if not np.isfinite(x) or x == 0:
    return x if np.isfinite(x) else None
  return float(f"{x:.{sig}g}")


def pct(a, ps=(5, 50, 95), mean=True, w=None):
  a = np.asarray(a, dtype=float)
  ok = np.isfinite(a)
  if w is not None:
    w = np.asarray(w, dtype=float)[ok]
  a = a[ok]
  if a.size == 0:
    return None
  if w is None:
    d = {f"p{p}": num(np.percentile(a, p)) for p in ps}
  else:                                   # weighted percentiles (length-weighted)
    o = np.argsort(a)
    c = np.cumsum(w[o]) / w.sum()
    d = {f"p{p}": num(a[o][np.searchsorted(c, p / 100)]) for p in ps}
  if mean:
    d["mean"] = num(np.average(a, weights=w))
  return d


def hist(a, edges):
  h, _ = np.histogram(a, bins=edges)
  return [int(x) for x in h]


def fiber_stats(X, E, P, L, fe, fid):
  """Per-fibre quantities from the fibre edges (chains)."""
  nf = fid.max() + 1
  arc = np.bincount(fid, weights=L[fe], minlength=nf)
  n_seg = np.bincount(fid, minlength=nf)

  # chain ends: nodes of fibre-degree 1, two per fibre for a simple chain
  Ef = E[fe]
  deg = np.bincount(Ef.ravel(), minlength=len(X))
  node_f = np.full(len(X), -1)
  node_f[Ef[:, 0]] = fid
  node_f[Ef[:, 1]] = fid
  ends = np.flatnonzero(deg == 1)
  ends = ends[np.argsort(node_f[ends], kind="stable")]
  ef = node_f[ends]
  n_ends = np.bincount(ef, minlength=nf)
  simple = n_ends == 2
  first = np.searchsorted(ef, np.arange(nf))
  if len(ends) == 0:   # no chain ends at all (synthetic grid: every node is a junction)
    a = b = np.zeros(nf, dtype=int)
  else:
    a = ends[np.minimum(first, len(ends) - 1)]
    b = ends[np.minimum(first + 1, len(ends) - 1)]
  chord = X[b] - X[a]
  proj = np.linalg.norm(chord, axis=1)
  shape = np.where(simple, proj / np.maximum(arc, 1e-300), np.nan)

  # in-plane orientation gamma in (-pi/2, pi/2], fibre chord (x = MD)
  gamma = np.arctan2(chord[:, 1], chord[:, 0])
  gamma = (gamma + np.pi / 2) % np.pi - np.pi / 2
  gamma = np.where(simple, gamma, np.nan)

  # segment orientation (length-weighted), and out-of-plane inclination
  t = X[Ef[:, 1]] - X[Ef[:, 0]]
  g_seg = np.arctan2(t[:, 1], t[:, 0])
  g_seg = (g_seg + np.pi / 2) % np.pi - np.pi / 2
  incl = np.degrees(np.arcsin(np.clip(np.abs(t[:, 2]) / np.maximum(L[fe], 1e-300), 0, 1)))

  # z-extent of each fibre
  zmin = np.full(nf, np.inf)
  zmax = np.full(nf, -np.inf)
  for c in (0, 1):
    np.minimum.at(zmin, fid, X[Ef[:, c], 2])
    np.maximum.at(zmax, fid, X[Ef[:, c], 2])

  # eta of f(gamma) = (1 + eta cos 2 gamma) / pi:  E[cos 2 gamma] = eta / 2
  ok = np.isfinite(gamma)
  eta_chord = 2 * np.mean(np.cos(2 * gamma[ok]))
  eta_seg = 2 * np.average(np.cos(2 * g_seg), weights=L[fe])

  # cross-section columns: one value per fibre (first edge), coarseness = mass / length
  Pf = P[fe]
  coarse = np.bincount(fid, weights=Pf[:, 0], minlength=nf) / np.maximum(arc, 1e-300)
  first_edge = np.unique(fid, return_index=True)[1]
  cols = {}
  for c, name in [(1, "EA"), (2, "kG_1A"), (4, "G_xI_x"), (5, "E_1I_1"), (6, "E_2I_2"),
                  (13, "width1"), (14, "width2")]:
    cols[name] = pct(Pf[first_edge, c], mean=False)
  cols["kG_1A/EA"] = pct(Pf[first_edge, 2] / Pf[first_edge, 1], mean=False)
  # per-fibre constancy: max relative spread of EA along a fibre
  ea_max = np.zeros(nf)
  ea_min = np.full(nf, np.inf)
  np.maximum.at(ea_max, fid, Pf[:, 1])
  np.minimum.at(ea_min, fid, Pf[:, 1])
  cols["EA_spread_within_fiber_max"] = num(np.nanmax(ea_max / ea_min - 1))
  cols["distinct_EA_values"] = int(len(np.unique(Pf[:, 1])))
  cols["coarseness_ug_per_m"] = pct(coarse * 1e15)      # kg/um -> ug/m

  return {
    "count": int(nf),
    "simple_chains": int(simple.sum()),
    "arc_length": pct(arc),
    "projected_length": pct(proj[simple]),
    "shape_factor": pct(shape),
    "segments_per_fiber": pct(n_seg),
    "segment_length": pct(L[fe]),
    "orientation": {
      "eta_cos1_from_chords": num(eta_chord),
      "eta_cos1_from_segments_lenweighted": num(eta_seg),
      "gamma_deg_hist_-90_90_by_15": hist(np.degrees(gamma[ok]), np.arange(-90, 91, 15)),
    },
    "inclination_deg_lenweighted": pct(incl, ps=(50, 90, 99), w=L[fe]),
    "z_extent": pct(zmax - zmin),
    "columns": cols,
  }, arc


def bond_stats(X, E, P, L, be, fe, fid_edge, node_f, arc):
  Eb = E[be]
  fa, fb = node_f[Eb[:, 0]], node_f[Eb[:, 1]]
  d = X[Eb[:, 1]] - X[Eb[:, 0]]
  dz = np.abs(d[:, 2])
  dxy = np.linalg.norm(d[:, :2], axis=1)
  Lb = L[be]
  incl = np.degrees(np.arctan2(dz, dxy))

  pair = np.sort(np.c_[fa, fb], axis=1)
  key = pair[:, 0].astype(np.int64) * (node_f.max() + 2) + pair[:, 1]
  ukey, first, inv, cnt = np.unique(key, return_index=True, return_inverse=True, return_counts=True)

  # nearest other bond of the same fibre pair (midpoint distance)
  mid = 0.5 * (X[Eb[:, 0]] + X[Eb[:, 1]])
  o = np.argsort(inv, kind="stable")
  bounds = np.r_[0, np.cumsum(cnt)]
  nn = np.full(len(Eb), np.nan)
  for g in np.flatnonzero(cnt > 1):
    idx = o[bounds[g]:bounds[g + 1]]
    m = mid[idx]
    D = np.linalg.norm(m[:, None] - m[None], axis=2)
    np.fill_diagonal(D, np.inf)
    nn[idx] = D.min(1)

  bonds_per_node = np.bincount(Eb.ravel(), minlength=len(X))
  on_fiber = node_f >= 0
  nf = len(arc)
  bonds_per_fiber = np.bincount(np.r_[fa, fb], minlength=nf)
  upair = pair[first]
  partners = np.bincount(upair.ravel(), minlength=nf)
  Pb = P[be]
  return {
    "count": int(be.sum()),
    "same_fiber": int((fa == fb).sum()),
    "end_not_on_fiber": int(((fa < 0) | (fb < 0)).sum()),
    "length": pct(Lb),
    "abs_dz": pct(dz),
    "dxy": pct(dxy),
    "inclination_deg": pct(incl),
    "frac_steeper_than_45deg": num((incl > 45).mean()),
    "bonded_fiber_pairs": int(len(ukey)),
    "bonds_per_pair_hist_1_to_8plus": [int((cnt == k).sum()) for k in range(1, 8)] + [int((cnt >= 8).sum())],
    "bonds_per_pair": pct(cnt, ps=(50, 90, 99)),
    "spacing_same_pair_nn": pct(nn, ps=(5, 25, 50, 75, 95)),
    "per_fiber": pct(bonds_per_fiber),
    "per_mm_fiber": num(be.sum() * 2 / arc.sum() * 1e3),
    "partners_per_fiber": pct(partners),
    "fibers_without_bond": int((bonds_per_fiber == 0).sum()),
    "per_node": {"max": int(bonds_per_node.max()),
                 "frac_fiber_nodes_with_bond": num((bonds_per_node[on_fiber] > 0).mean()),
                 "hist_0_to_5plus": [int((bonds_per_node[on_fiber] == k).sum()) for k in range(5)]
                                    + [int((bonds_per_node[on_fiber] >= 5).sum())]},
    "columns": {name: pct(Pb[:, c], mean=False)
                for c, name in [(0, "mass"), (1, "EA"), (2, "kG_1A"), (4, "G_xI_x"),
                                (5, "E_1I_1"), (6, "E_2I_2"), (13, "width1"), (14, "width2")]},
  }


def z_stats(X, E, L, be, fe, lo, hi, nbins=10):
  T = hi[2] - lo[2]
  zr = (X[:, 2] - lo[2]) / T
  mz = 0.5 * (zr[E[:, 0]] + zr[E[:, 1]])
  edges = np.linspace(0, 1, nbins + 1)
  fl, _ = np.histogram(mz[fe], bins=edges, weights=L[fe])
  bl, _ = np.histogram(mz[be], bins=edges)

  def cut(axis, frac=0.5):
    c = lo[axis] + frac * (hi[axis] - lo[axis])
    s = (X[E[:, 0], axis] - c) * (X[E[:, 1], axis] - c) < 0
    other = [i for i in range(3) if i != axis]
    area = np.prod([hi[i] - lo[i] for i in other])
    return s, area

  out = {"fiber_length_frac_per_tenth": [num(x) for x in fl / fl.sum()],
         "bonds_frac_per_tenth": [num(x) for x in bl / max(bl.sum(), 1)]}
  cuts = {}
  for name, ax in (("x", 0), ("y", 1), ("z", 2)):
    s, area = cut(ax)
    cuts[name] = {"fiber_per_mm2": num((s & fe).sum() / area * 1e6),
                  "bond_per_mm2": num((s & be).sum() / area * 1e6)}
  out["midplane_cuts"] = cuts
  zc = []
  for frac in (0.25, 0.5, 0.75):
    s, area = cut(2, frac)
    zc.append(num(s.sum() / area * 1e6))
  out["z_cuts_per_mm2_at_quarter_half_threequarter"] = zc
  cx = (cuts["x"]["fiber_per_mm2"] or 0) + (cuts["x"]["bond_per_mm2"] or 0)
  cz = (cuts["z"]["fiber_per_mm2"] or 0) + (cuts["z"]["bond_per_mm2"] or 0)
  out["anisotropy_cut_z_over_cut_x"] = num(cz / cx) if cx else None
  return out


def metis_stats(X, E, ks):
  """k-way METIS partitions: edge cut and per-axis variance explained by the parts."""
  exe = os.environ.get("GPMETIS") or shutil.which("gpmetis")
  if exe is None:
    tprint("gpmetis not found (set $GPMETIS), skipping partition")
    return None
  n = len(X)
  e = np.sort(E, axis=1)
  e = np.unique(e[e[:, 0] != e[:, 1]], axis=0)
  A = sp.coo_matrix((np.ones(2 * len(e), np.int8), (np.r_[e[:, 0], e[:, 1]], np.r_[e[:, 1], e[:, 0]])),
                    shape=(n, n)).tocsr()
  out = {}
  with tempfile.TemporaryDirectory() as tmp:
    gfile = os.path.join(tmp, "graph")
    with open(gfile, "w") as f:
      f.write(f"{n} {len(e)}\n")
      for r in range(n):
        f.write(" ".join(map(str, A.indices[A.indptr[r]:A.indptr[r + 1]] + 1)) + "\n")
    for k in ks:
      tprint(f"gpmetis k = {k}")
      subprocess.run([exe, gfile, str(k)], check=True, capture_output=True)
      part = np.loadtxt(f"{gfile}.part.{k}", dtype=np.int64)
      cut = int((part[e[:, 0]] != part[e[:, 1]]).sum())
      cnt = np.bincount(part, minlength=k)
      r2 = {}
      for ax, name in enumerate("xyz"):
        x = X[:, ax]
        means = np.bincount(part, weights=x, minlength=k) / np.maximum(cnt, 1)
        tot = ((x - x.mean()) ** 2).sum()
        r2[name] = num((cnt * (means - x.mean()) ** 2).sum() / tot) if tot > 0 else None
      out[f"k{k}"] = {"edge_cut": cut, "r2": r2,
                      "part_z_mean": [num(m) for m in
                                      np.bincount(part, weights=X[:, 2], minlength=k) / np.maximum(cnt, 1)]}
  return out


def net_stats(path, metis=()):
  tprint(f"reading {path}")
  with h5py.File(path, "r") as f:
    X = f["domain/points"][:]
    E = f["domain/edges"][:].astype(np.int64)
    P = f["domain/properties"][:]
    attrs = {k: (v.tolist() if hasattr(v, "tolist") else v) for k, v in f["domain"].attrs.items()}

  L = np.linalg.norm(X[E[:, 1]] - X[E[:, 0]], axis=1)
  fid_edge = P[:, 15].astype(np.int64)
  be = fid_edge == -1
  fe = ~be
  # compact fibre ids 0..nf-1 (pruning in make_geo2 can leave gaps)
  _, fid = np.unique(fid_edge[fe], return_inverse=True)
  fid_c = np.full(len(E), -1)
  fid_c[fe] = fid
  node_f = np.full(len(X), -1)
  node_f[E[fe, 0]] = fid
  node_f[E[fe, 1]] = fid

  lo, hi = X.min(0), X.max(0)
  ext = hi - lo
  mass = np.nansum(P[:, 0]) if np.isfinite(P[fe, 0]).any() else np.nan
  area = ext[0] * ext[1]
  gsm = mass / area * 1e15                   # kg/um^2 -> g/m^2
  sheet = {
    "nodes": int(len(X)), "edges": int(len(E)),
    "fiber_edges": int(fe.sum()), "bond_edges": int(be.sum()),
    "extent_um": [num(x) for x in ext],
    "grammage_gsm": num(gsm),
    "density_kg_m3": num(gsm * 1e-3 / (ext[2] * 1e-6)),
    "fiber_length_per_area_um_per_um2": num(L[fe].sum() / area),
    "info": {k: v for k, v in attrs.items() if k not in ("units",)},
  }
  tprint("fibres")
  fs, arc = fiber_stats(X, E, P, L, fe, fid)
  tprint("bonds")
  bs = bond_stats(X, E, P, L, be, fe, fid_c, node_f, arc)
  tprint("z")
  zs = z_stats(X, E, L, be, fe, lo, hi)
  doc = {"sheet": sheet, "fibers": fs, "bonds": bs, "z": zs}
  if metis:
    doc["partition"] = metis_stats(X, E, metis)
  return doc


def main():
  ap = argparse.ArgumentParser(description="net_stats by Joseph Holten")
  ap.add_argument("inputs", nargs="+", help=".geo.h5 domains (make_geo2 format)")
  ap.add_argument("--metis", type=int, nargs="+", default=[], metavar="K",
                  help="also report METIS k-way partitions for these K")
  args = ap.parse_args()
  for path in args.inputs:
    doc = {os.path.basename(path).removesuffix(".geo.h5"): net_stats(path, args.metis)}
    yaml.safe_dump(doc, sys.stdout, sort_keys=False, default_flow_style=None, width=120)
    sys.stdout.flush()


if __name__ == "__main__":
  main()
