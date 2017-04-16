import numpy as np
import numpy.linalg

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

def gramschmidtw(s, s0=(), E=()): # scalar product L2 space
  if s0 == ():
    s0 = np.ones(s.n) # constant function
  n = len(s0)
  if E == ():
    E = np.diagflat(np.ones(n)) # identity
  S = np.empty((n, n))
  S[0] = s0 / np.sqrt(sum(s0**2 * s.w))
  for k in range(1, n):
    S[k] = E[k]
    for j in range(k):
      S[k] = S[k] - sum(E[k] * S[j] * s.w) * S[j]
    S[k] = S[k] / np.sqrt(sum(S[k]**2 * s.w))
  return S.T
gramschmidtBoundary = gramschmidtw

def gramschmidt(s0, E=()): # needs n
  n = len(s0)
  if E == ():
    E = np.diagflat(np.ones(n))
  S = np.empty((n, n))
  S[0] = s0 / numpy.linalg.norm(s0)
  for k in range(1, n):
    S[k] = E[k]
    for j in range(k):
      S[k] = S[k] - E[k].dot(S[j]) * S[j]
    S[k] = S[k] / numpy.linalg.norm(S[k])
  return S.T
