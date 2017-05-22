from __load__ import *
a = 3.
b = 3.
nso = 10


def ex_k(t, s, a, b):
  return - a*b / np.pi / (a**2 + b**2 - (a**2 - b**2) * np.cos( t + s )) * (2 * np.pi / len(s))

def ex_w(n):
  d = np.concatenate(([4], 4 * np.ones(n - 2) ,[4]))
  m = (np.diagflat(d) + np.diagflat(np.ones(n-1), 1) + np.diagflat(np.ones(n-1), -1) )
  m[0, -1] = 1
  m[-1, 0] = 1
  m = m / 6
  return m

def ex_phi(t):
  # return np.exp(np.cos(t)) * np.cos(np.sin(t))
  return t**3
def ex_f(t, a, b):
  c = (a - b) / (a + b)
  # return ex_phi(t) + np.exp(c * np.cos(t)) * np.cos(c * np.sin(t))
  return ex_phi(t) + (2 * np.pi)**4 / 4 * (- a*b / np.pi / (a**2 + b**2 - (a**2 - b**2) * np.cos( 0 )))


def func(nso):
  so = sg.Segment(nso, f_inargs = (sh.ellipse, (0, a, b)), quad='ps')
  Kp = ly.layerpotSD(s=so)
  K = ly.layerpotD(s=so)

  s = 2 * np.pi * np.array([so.t for k in range(len(so.t))])
  t = s.T

  phi = ex_phi(2 * np.pi * so.t)
  com_f = -ex_k(t, s, a, b).dot(phi) + phi
  f = ex_f(2 * np.pi * so.t, a, b)
  print(nso)
  return max(abs(com_f - f))
def func2(nso):
  so = sg.Segment(nso, f_inargs = (sh.ellipse, (0, a, b)), quad='ps')
  Kp = ly.layerpotSD(s=so)
  K = ly.layerpotD(s=so)

  s = 2 * np.pi * np.array([so.t for k in range(len(so.t))])
  t = s.T

  q = so.t * 2 * np.pi
  phi = ex_phi(2 * np.pi * so.t)
  com_f = -ex_k(t, s, a, b).dot(phi) + phi
  f = ex_f(2 * np.pi * so.t, a, b)
  print(nso)
  return (q, phi, k, com_f)
def func3(nso):
  so = sg.Segment(nso, f_inargs = (sh.ellipse, (0, a, b)), quad='ps')
  Kp = ly.layerpotSD(s=so)
  K = ly.layerpotD(s=so)

  s = 2 * np.pi * so.t
  t = 0

  q = 2 * np.pi * so.t
  phi = ex_phi(2 * np.pi * so.t)
  k = -ex_k(t, s, a, b)
  com_f = k.dot(phi) + ex_phi(t * 2 * np.pi)
  f = ex_f(2 * np.pi * t, a, b)
  print(nso)
  # return (q, phi, k, com_f, t)
  return abs(com_f - f)
def func_it(ntot):
  err = np.ones(ntot, float)
  n = np.array([4*2**k for k in range(ntot)])
  n = np.array([4 + k for k in range(ntot)])
  for k in range(ntot):
    err[k] = func(n[k])
  return err
