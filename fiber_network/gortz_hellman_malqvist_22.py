import numpy as np
import scipy.sparse as sp

from joblib import Parallel, delayed


def _coarse_basis_2d(points, n_elem_1d, epsilon=1e-10):
  # B contains (N+1)^2 columns and length(p) rows.
  # Each column represents a bilinear basis function on the mesh of N*N squares of side length H=1/N

  min_coord = np.min(points, axis=0)
  h = np.divide( (np.max(points, axis=0) - min_coord)[:len(n_elem_1d)], n_elem_1d )

  m, n = (n_elem_1d[0]+1) * (n_elem_1d[1]+1), len(points)  # Dimensions of the resulting matrix
  node_vec, index_vec, value_vec = [], [], []

  # loop over elements to find nodes
  for i in range(n_elem_1d[0]):
    for j in range(n_elem_1d[1]):
      nodes = np.where( (points[:,0] - min_coord[0] + epsilon > h[0] * i) &
                        (points[:,0] - min_coord[0] + epsilon < h[0] * (i + 1)) &
                        (points[:,1] - min_coord[1] + epsilon > h[1] * j) &
                        (points[:,1] - min_coord[1] + epsilon < h[1] * (j + 1)) )[0]
      index = [ j*(n_elem_1d[0]+1)+i,       j*(n_elem_1d[0]+1)+i+1,
                (j+1)*(n_elem_1d[0]+1)+i+1, (j+1)*(n_elem_1d[0]+1)+i ]
      x = (points[nodes,0] - min_coord[0]) / h[0] - i
      y = (points[nodes,1] - min_coord[1]) / h[1] - j

      node_vec.extend( [nodes, nodes, nodes, nodes] )
      index_vec.extend( [index[0]*np.ones(len(nodes)), index[1]*np.ones(len(nodes)),
                         index[2]*np.ones(len(nodes)), index[3]*np.ones(len(nodes))] )
      value_vec.extend( [(1-x)*(1-y), x*(1-y), x*y, (1-x)*y] )

  node_vec, index_vec, value_vec = np.hstack(node_vec), np.hstack(index_vec), np.hstack(value_vec)
  return sp.csr_matrix((value_vec, (node_vec, index_vec)), shape=(n, m))


def _coarse_basis_3d(points, n_elem_1d, epsilon=1e-10):
  # B contains (N+1)^3 columns and length(p) rows.
  # Each column represents a bilinear basis function on the mesh of N*N squares of side length H=1/N

  min_coord = np.min(points, axis=0)
  h = np.divide( (np.max(points, axis=0) - min_coord)[:len(n_elem_1d)], n_elem_1d )

  m, n = (n_elem_1d[0]+1) * (n_elem_1d[1]+1) * (n_elem_1d[2]+1), len(points)
  node_vec, index_vec, value_vec = [], [], []

  # loop over elements to find nodes
  for i in range(n_elem_1d[0]):
    for j in range(n_elem_1d[1]):
      for k in range(n_elem_1d[2]):
        nodes = np.where( (points[:,0] - min_coord[0] + epsilon > h[0] * i) &
                          (points[:,0] - min_coord[0] + epsilon < h[0] * (i + 1)) &
                          (points[:,1] - min_coord[1] + epsilon > h[1] * j) &
                          (points[:,1] - min_coord[1] + epsilon < h[1] * (j + 1)) &
                          (points[:,2] - min_coord[2] + epsilon > h[2] * k) &
                          (points[:,2] - min_coord[2] + epsilon < h[2] * (k + 1)) )[0]
        index = [ (k * (n_elem_1d[1]+1) + j) * (n_elem_1d[0]+1) + i,
                  (k * (n_elem_1d[1]+1) + j) * (n_elem_1d[0]+1) + i+1,
                  (k * (n_elem_1d[1]+1) + j+1) * (n_elem_1d[0]+1) + i+1,
                  (k * (n_elem_1d[1]+1) + j+1) * (n_elem_1d[0]+1) + i,
                  ((k+1) * (n_elem_1d[1]+1) + j) * (n_elem_1d[0]+1) + i,
                  ((k+1) * (n_elem_1d[1]+1) + j) * (n_elem_1d[0]+1) + i+1, 
                  ((k+1) * (n_elem_1d[1]+1) + j+1) * (n_elem_1d[0]+1) + i+1,
                  ((k+1) * (n_elem_1d[1]+1) + j+1) * (n_elem_1d[0]+1) + i ]
        x = (points[nodes,0] - min_coord[0]) / h[0] - i
        y = (points[nodes,1] - min_coord[1]) / h[1] - j
        z = (points[nodes,2] - min_coord[2]) / h[2] - k

        node_vec.extend( [nodes, nodes, nodes, nodes, nodes, nodes, nodes, nodes] )
        index_vec.extend( [index[0]*np.ones(len(nodes)), index[1]*np.ones(len(nodes)),
                           index[2]*np.ones(len(nodes)), index[3]*np.ones(len(nodes)),
                           index[4]*np.ones(len(nodes)), index[5]*np.ones(len(nodes)),
                           index[6]*np.ones(len(nodes)), index[7]*np.ones(len(nodes))] )
        value_vec.extend( [ (1-z)*(1-x)*(1-y), (1-z)*x*(1-y), (1-z)*x*y, (1-z)*(1-x)*y,
                            z*(1-x)*(1-y), z*x*(1-y), z*x*y, z*(1-x)*y ] )

  node_vec, index_vec, value_vec = np.hstack(node_vec), np.hstack(index_vec), np.hstack(value_vec)
  return sp.csr_matrix((value_vec, (node_vec, index_vec)), shape=(n, m))


class gortz_hellman_malqvist_22:
  def __init__( self, points, n_elem_1d, epsilon=1e-10, repeat=1 ):
    coarse_basis, int_nodes = [], []
    if   len(n_elem_1d) == 2:
      coarse_basis = _coarse_basis_2d(points, n_elem_1d,  epsilon)
      int_nodes    = np.concatenate([ j*(n_elem_1d[0]+1) + np.arange(1, n_elem_1d[0])
                                      for j in range(1, n_elem_1d[1]) ])  # find interior nodes
    elif len(n_elem_1d) == 3:
      coarse_basis = _coarse_basis_3d(points, n_elem_1d, epsilon)
      int_nodes    = np.concatenate([
        (k * (n_elem_1d[1]+1)+j) * (n_elem_1d[0]+1) + np.arange(1, n_elem_1d[0])
        for j in range(1, n_elem_1d[1]) for k in range(1, n_elem_1d[2]) ])  # find interior nodes

    coarse_basis_int = coarse_basis[:, int_nodes]
    bnd_nodes = np.where(np.sum(coarse_basis_int, axis=1) < epsilon)[0]  # coarse basis in V

    # Get the row indices and data of the column to be set to zero
    for index in bnd_nodes:
      start_idx, end_idx = coarse_basis.indptr[index], coarse_basis.indptr[index + 1]
      coarse_basis.data[start_idx:end_idx] = 0.

    if repeat > 1:
      helper = np.eye(repeat)
      coarse_basis = sp.kron(coarse_basis, helper)
      coarse_basis_int = sp.kron(coarse_basis_int, helper)

    self.coarse_basis     = sp.csr_matrix(coarse_basis)
    self.coarse_basis_int = sp.csr_matrix(coarse_basis_int)

    # cached factorizations, built on the first precond call and reused across all
    # CG iterations (see _ensure_setup)
    self._solver_id = None  # id() of the lhs_mat the cache was built for
    self._supports  = None  # per-subdomain support row indices


  def _ensure_setup(self, lhs_mat, epsilon):
    # Factorize the coarse and local operators once and reuse them for every
    # precond call that shares the same lhs_mat (all CG iterations). This is the
    # same LU factorization spsolve does internally -- only the per-call
    # refactorization is removed, so the preconditioner applied is unchanged.
    if self._solver_id == id(lhs_mat):
      return

    # support (rows with weight > epsilon) of each coarse basis column defines one
    # subdomain; depends only on coarse_basis, so compute it once.
    if self._supports is None:
      cb = self.coarse_basis.tocsc()
      self._supports = [ np.sort(cb.indices[cb.indptr[k]:cb.indptr[k+1]]
                                 [cb.data[cb.indptr[k]:cb.indptr[k+1]] > epsilon])
                         for k in range(cb.shape[1]) ]

    coarse_op = self.coarse_basis_int.T.dot(lhs_mat.dot(self.coarse_basis_int))
    self._coarse_solve = sp.linalg.factorized(sp.csc_matrix(coarse_op))
    self._local_solve  = [ sp.linalg.factorized(sp.csc_matrix(lhs_mat[nj, :][:, nj]))
                           if nj.size else None
                           for nj in self._supports ]
    self._solver_id = id(lhs_mat)


  def precond(self, lhs_mat, rhs_vec, n_jobs=None, epsilon=1e-14):
    self._ensure_setup(lhs_mat, epsilon)

    # coarse correction: B_int (B_int^T A B_int)^-1 B_int^T rhs
    result_vec = self.coarse_basis_int.dot(
      self._coarse_solve(self.coarse_basis_int.T.dot(rhs_vec)))

    # local corrections: the restriction Ij^T rhs is just rhs[nj] and the
    # prolongation Ij z is a scatter-add result[nj] += z, so the dense n x |nj|
    # matrix Ij is never formed.
    if n_jobs is None:
      for nj, solve in zip(self._supports, self._local_solve):
        if solve is not None:
          result_vec[nj] += solve(rhs_vec[nj])
      return result_vec

    # threads, not processes: the cached factorizations are not picklable, and the
    # scatter-add is kept serial below to avoid races on result_vec.
    def process(nj, solve):
      return nj, solve(rhs_vec[nj])

    result = Parallel(n_jobs=n_jobs, prefer="threads")(
      delayed(process)(nj, solve)
      for nj, solve in zip(self._supports, self._local_solve) if solve is not None)
    for nj, z in result:  result_vec[nj] += z

    return result_vec
