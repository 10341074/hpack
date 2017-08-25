import numpy as np
import numpy.linalg

def trigonometricBasis(n):
  B = np.empty((n, n), float)
  t = 2 * np.pi * np.arange(n) / n
  if n % 2 == 1:
    print('WARNING: not even trigonometric basis')
  B[:, 0] = 1
  for j in range(1, n):
    k = np.ceil(j/2)
    if  j % 2 == 1:
      B[:, j] = np.cos(k * t)
    else:
      B[:, j] = np.sin(k * t)
  return B

def eigmaxpower(A, nit=100):
  it = 0
  q = np.ones(len(A))
  while it < nit:
    z = A.dot(q)
    q = z / numpy.linalg.norm(z)
    nu = q.T.dot(A.dot(q))
    it = it + 1
  return (nu, q)

# check it's the same: q/qw constant, nu/nuw
def eigmaxpowerw(A, s, nit=100):
  it = 0
  q = np.ones(len(A))
  while it < nit:
    z = A.dot(q)
    q = z / np.sqrt(sum((z**2 * s.w)))
    nu = sum(q * A.dot(q) * s.w)
    it = it + 1
  return (nu, q)

def gramschmidtw(s=(), s0=(), E=(), w=()): # scalar product L2 space
  if w == ():
    w = s.w
  if s0 == ():
    s0 = np.ones(s.n) # constant function
  n = len(s0)
  if E == ():
    E = np.diagflat(np.ones(n)) # identity
  S = np.empty((n, n))
  S[0] = s0 / np.sqrt(sum(s0**2 * w))
  for k in range(1, n):
    S[k] = E[k]
    for j in range(k):
      S[k] = S[k] - sum(E[k] * S[j] * w) * S[j]
    S[k] = S[k] / np.sqrt(sum(S[k]**2 * w))
  return S.T
gramschmidtBoundary = gramschmidtw

def gramschmidt(s0, E=()): # needs n
  n = len(s0)
  if E == ():
    E = np.eye(n)
  S = np.empty((n, n))
  S[0] = s0 / numpy.linalg.norm(s0)
  for k in range(1, n):
    S[k] = E[k]
    for j in range(k):
      S[k] = S[k] - E[k].dot(S[j]) * S[j]
    S[k] = S[k] / numpy.linalg.norm(S[k])
  return S.T

def base_norm(A, w=()):
  n = A.shape[1]
  if w == ():
    w = np.ones(A.shape[0])
  for k in range(n):
    norm = np.sqrt(sum(A[:, k]**2 * w))
    A[:, k] = A[:, k] / norm
  return A

def base_mean(A, w=()):
  n = A.shape[1]
  if w == ():
    w = np.ones(A.shape[0])
  for k in range(n):
    mean = sum(A[:, k] * w) / sum(w)
    A[:, k] = A[:, k] - mean
  return A

# test one type of basis
def m_basis(so):
   so.BX = np.eye(so.n)
   so.BX = base_mean(so.BX, so.w)
   return so.BX
def m_change_basis(v):
  n = so.n
  M = np.eye(so.n)
  m = so.w.dot( np.eye(so.n) ) / sum(so.w)
  S = base_mean(M, so.w)

  return A
