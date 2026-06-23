#!/usr/bin/env python

"""
partition_net.py — partition a fiber-network .geo.h5 with METIS and store the
partition back into the file so net_vis.py can color by it.

The graph mirrors the one HyperHDG partitions in distribute_domain.hxx for
hyEdge_dim == 1 networks: vertices are the graph nodes (== hypernodes ==
points), graph edges are the beams (hyperedges), and each vertex is weighted by
its incident-edge count so every part gets a balanced share of assembly work.
We use serial METIS (METIS_PartGraphKway) here — for an offline visualization
step the parallel ParMETIS buys nothing over serial METIS, and the partition
quality is the same.

Two datasets are written under /domain and wired into the existing /VTKHDF view:
  node_partition  (PointData, per point) — the part each node was assigned to
  edge_partition  (CellData,  per beam, exposed as 'partition') — the beam's
                  owning part = min(part[n0], part[n1]), matching the
                  owner = min part rule in distribute_domain.hxx

Visualize with:
    net_vis.py net.geo.h5 --color-by partition --warp-by none --arrows 0

Inputs
------
A .geo.h5 written by make_geo.py: /domain/{points,edges} plus a /VTKHDF view.

Author
------
Joseph Holten, KIT, 2026.
"""

import argparse
import ctypes
import os
import sys
import time
from ctypes.util import find_library

import h5py
import numpy as np
import scipy.sparse as sp


def tprint(*args, **kwargs):
    print(f"[{time.strftime('%H:%M:%S')}]", *args, **kwargs)


def build_graph(edges, n_nodes):
    """Undirected simple-graph CSR (no self-loops, parallel edges collapsed).

    Returns (xadj, adjncy, vwgt) with vwgt = per-node incident (real) edge
    count, counted with multiplicity to match distribute_domain.hxx's degree
    weighting; isolated nodes are clamped to weight 1 so METIS stays happy.
    """
    i, j = edges[:, 0], edges[:, 1]
    keep = i != j  # drop self-loops from the adjacency
    ii, jj = i[keep], j[keep]

    # symmetric adjacency; sum_duplicates then flatten to unweighted presence
    rows = np.concatenate([ii, jj])
    cols = np.concatenate([jj, ii])
    A = sp.coo_matrix((np.ones(rows.shape[0], dtype=np.int8), (rows, cols)),
                      shape=(n_nodes, n_nodes)).tocsr()
    A.sum_duplicates()
    A.data[:] = 1

    xadj = A.indptr.astype(np.int32)
    adjncy = A.indices.astype(np.int32)
    # vwgt: incident-edge count with multiplicity (work ~ real beams touched)
    vwgt = np.bincount(rows, minlength=n_nodes).astype(np.int32)
    vwgt = np.maximum(vwgt, 1)
    return xadj, adjncy, vwgt


def find_metis_lib(explicit=None):
    """Locate libmetis.so: explicit path, $METIS_LIB, or the loader's search."""
    for cand in (explicit, os.environ.get("METIS_LIB"), find_library("metis")):
        if cand and os.path.exists(cand):
            return cand
        if cand and not os.path.sep in cand:
            return cand  # let the dynamic loader resolve a bare soname
    return None


def partition_metis_ctypes(libpath, xadj, adjncy, vwgt, nparts):
    """Call METIS_PartGraphKway via ctypes. METIS here is 32-bit idx/real."""
    lib = ctypes.CDLL(libpath)
    idx_t = ctypes.c_int32   # IDXTYPEWIDTH == 32
    real_t = ctypes.c_float  # REALTYPEWIDTH == 32
    IDXP, REALP = ctypes.POINTER(idx_t), ctypes.POINTER(real_t)

    xadj = np.ascontiguousarray(xadj, dtype=np.int32)
    adjncy = np.ascontiguousarray(adjncy, dtype=np.int32)
    vwgt = np.ascontiguousarray(vwgt, dtype=np.int32)
    n = xadj.shape[0] - 1
    part = np.empty(n, dtype=np.int32)

    nvtxs, ncon = idx_t(n), idx_t(1)
    nparts_c, objval = idx_t(nparts), idx_t(0)

    lib.METIS_PartGraphKway.restype = ctypes.c_int
    lib.METIS_PartGraphKway.argtypes = [IDXP] * 4 + [IDXP] * 3 + \
        [IDXP, REALP, REALP, IDXP, IDXP, IDXP]

    def p(a):
        return a.ctypes.data_as(IDXP)

    # options = NULL -> METIS defaults, incl. 0-based (C) numbering
    ret = lib.METIS_PartGraphKway(
        ctypes.byref(nvtxs), ctypes.byref(ncon), p(xadj), p(adjncy),
        p(vwgt), None, None,
        ctypes.byref(nparts_c), None, None, None,
        ctypes.byref(objval), p(part))
    if ret != 1:  # METIS_OK == 1
        raise RuntimeError(f"METIS_PartGraphKway failed (code {ret})")
    return part, int(objval.value)


def partition_pymetis(xadj, adjncy, vwgt, nparts):
    import pymetis
    edgecut, membership = pymetis.part_graph(
        nparts, xadj=xadj.tolist(), adjncy=adjncy.tolist(),
        vweights=vwgt.tolist())
    return np.asarray(membership, dtype=np.int32), int(edgecut)


def partition(xadj, adjncy, vwgt, nparts, metis_lib=None):
    """Partition into nparts; nparts==1 is trivial. Prefer libmetis, else pymetis."""
    n = xadj.shape[0] - 1
    if nparts <= 1:
        return np.zeros(n, dtype=np.int32), 0

    lib = find_metis_lib(metis_lib)
    if lib:
        tprint(f"partitioning with METIS_PartGraphKway via {lib}")
        return partition_metis_ctypes(lib, xadj, adjncy, vwgt, nparts)
    try:
        tprint("libmetis.so not found; falling back to pymetis")
        return partition_pymetis(xadj, adjncy, vwgt, nparts)
    except ImportError:
        sys.exit("error: no METIS available. Pass --metis-lib /path/to/libmetis.so, "
                 "set $METIS_LIB, or install pymetis.")


def quotient_coloring(node_part, edges, nparts):
    """Greedily color the quotient graph (parts as super-nodes, adjacent iff a
    beam crosses between them) so neighboring parts get different colors.

    Welsh-Powell: color parts in descending degree order, each taking the
    smallest color unused by its already-colored neighbors. Not guaranteed
    chromatic-minimal, but tight in practice (<=4 for planar-ish 2D networks).
    Returns (part_color, n_colors).
    """
    pa, pb = node_part[edges[:, 0]], node_part[edges[:, 1]]
    cut = pa != pb
    # unique unordered cross-part pairs -> tiny quotient graph
    pairs = np.unique(np.sort(np.stack([pa[cut], pb[cut]], axis=1), axis=1), axis=0)
    adj = [set() for _ in range(nparts)]
    for a, b in pairs:
        adj[a].add(int(b))
        adj[b].add(int(a))

    part_color = np.full(nparts, -1, dtype=np.int32)
    for p in sorted(range(nparts), key=lambda q: -len(adj[q])):
        used = {part_color[q] for q in adj[p] if part_color[q] >= 0}
        c = 0
        while c in used:
            c += 1
        part_color[p] = c
    n_colors = int(part_color.max()) + 1 if nparts else 0
    tprint(f"quotient graph: {nparts} parts, {len(pairs)} adjacent pairs, "
           f"colored with {n_colors} colors")
    return part_color, n_colors


def write_arrays(path, arrays, attrs=()):
    """Store int datasets under /domain and link them into the /VTKHDF view.

    arrays: list of (domain_name, data, assoc, link_name) where assoc is
    'PointData', 'CellData', or None (store only, no VTKHDF link).
    attrs: list of (name, value) written onto /domain.
    """
    tprint(f"writing {len(arrays)} datasets into '{path}'")
    with h5py.File(path, "a") as f:
        dom = f["domain"]
        if "VTKHDF" not in f:
            sys.exit("error: file has no /VTKHDF view; (re)run make_geo.py first")
        root = f["VTKHDF"]
        for dname, data, assoc, link in arrays:
            if dname in dom:
                del dom[dname]
            dom.create_dataset(dname, data=data.astype(np.int32), compression="gzip")
            if assoc:
                g = root.require_group(assoc)
                if link in g:
                    del g[link]
                g[link] = h5py.SoftLink(f"/domain/{dname}")
        for name, value in attrs:
            dom.attrs[name] = value


if __name__ == "__main__":
    p = argparse.ArgumentParser(description="partition a fiber-network .geo.h5 with METIS")
    p.add_argument("input", help="path to .geo.h5 (written by make_geo.py)")
    p.add_argument("-n", "--nparts", type=int, default=4, help="number of parts")
    p.add_argument("--metis-lib", default=None,
                   help="path to libmetis.so (else $METIS_LIB / loader search / pymetis)")
    args = p.parse_args()

    with h5py.File(args.input, "r") as f:
        points = f["domain/points"][...]
        edges = f["domain/edges"][...].astype(np.int64)
    n_nodes, n_edges = points.shape[0], edges.shape[0]
    tprint(f"read {n_nodes} nodes, {n_edges} edges from '{args.input}'")

    xadj, adjncy, vwgt = build_graph(edges, n_nodes)
    node_part, edgecut = partition(xadj, adjncy, vwgt, args.nparts,
                                   metis_lib=args.metis_lib)

    # beam owner = min part of its two endpoints (matches distribute_domain.hxx)
    edge_part = np.minimum(node_part[edges[:, 0]], node_part[edges[:, 1]])

    sizes = np.bincount(node_part, minlength=args.nparts)
    tprint(f"nparts={args.nparts}  edgecut={edgecut}")
    tprint(f"node counts per part: {sizes.tolist()}")
    tprint(f"  balance: max/avg = {sizes.max() / (n_nodes / args.nparts):.3f}")

    # quotient-graph coloring: few colors, adjacent parts always differ
    part_color, n_colors = quotient_coloring(node_part, edges, args.nparts)
    node_color = part_color[node_part]
    edge_color = part_color[edge_part]

    write_arrays(args.input, [
        ("node_partition", node_part,  "PointData", "node_partition"),
        ("edge_partition", edge_part,  "CellData",  "partition"),
        ("node_color",     node_color, "PointData", "node_color"),
        ("edge_color",     edge_color, "CellData",  "color"),
    ], attrs=[("n_parts", args.nparts), ("n_colors", n_colors)])
    tprint("done")
