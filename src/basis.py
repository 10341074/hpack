import numpy as np

def basis_lagr_trig(n = 10, t = ()):
  if t == ():
    t = np.linspace(0, 2*np.pi, 101)
  B = np.empty((2*n, len(t)))
  tn = np.pi / n * np.arange(2*n)
  for k, (tj, r) in enumerate(zip(tn, B)):
    B[k] = 1.0 / 2 / n * (1 + np.ones((1, n-2)).dot(np.cos(
        np.array(np.arange(1, n-1), ndmin=2).T.dot(
            np.array(t - tj, ndmin=2)
        )
    ))
    + np.cos(n * (t - tj))
    )
  return B.T

def basis_lagr(n = 10, t = ()):
  if t == ():
    t = np.linspace(0, 2*np.pi, 101)
  B = np.empty((2*n, len(t)))
  tn = np.pi / n * np.arange(2*n)
  B[0] = np.ones((1, len(t)))
  for k, tj in zip(np.arange(n), tn):
    B[2 * k - 1] = np.cos(k * (t - tj))
    B[2 * k] = np.sin(k * (t - tj))
  return B.T
