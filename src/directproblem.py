#!../../venv/bin/python
import matplotlib.pyplot as plt
import numpy as np
from scipy.linalg import solve

import my_type
import layerpot as ly
import segment as sg
import plot

def directA(l, c):
  Kp = l.exp(0, l.b, [], [])
  A = Kp
  A[np.diag_indices(len(l.b.x))]=A[np.diag_indices(len(l.b.x))] + 0.5 * c
  return A

def directrhs(b, z0):
  d = b.x - z0
  r = abs(d)
  cosphi = np.real(np.conj(b.nx) * d) / r
  rhs = - - ly.fundsol_deriv(r, cosphi, 0)
  return rhs
  
def directpb(z0, l, c, A=[], rhs=[]):
  if A == []:
    print 'Warning: calculated A, not passed'
    A = directA(l, c)
  if rhs == []:
    rhs = directrhs(l.b[0], z0)
  psi = solve(A, rhs)
  return psi

if __name__ == "__main__":
  h = 20
  n = 20
  c = 1. * (h + 1)/(h-1)

  # (K' + 0.5 * c(h) * I) psi = -phi_nu

  s = sg.Segment(n, f = sg.circle, inargs = (0, 1), periodic=True)
  # s = sg.Segment(n, Z=sg.kZ, Zp=sg.kZp, Zpp=sg.kZpp, args=[], periodic=True)
  b = sg.Boundary([s])
  p = sg.Pointset(s.x)

  l = sg.Layer(b, ly.layerpotSD)
  Kp = l.exp(0, l.b, [], [])
  A = Kp
  A[np.diag_indices(len(s.x))]=A[np.diag_indices(len(s.x))] + 0.5 * c
  A = directA(l, c)
  z0 = 2.1
  d = s.x - z0
  r = abs(d)
  cosphi = np.real(np.conj(s.nx) * d) / r
  # rhs = - ly.fundsol_deriv(r, cosphi, 0)
  rhs = directrhs(s, z0)
  psi = solve(A, rhs)


  x, y, pp = plot.meshgrid((-2, 2, 80))
  G = sg.Layer(b=b, exp=ly.layerpotS, dns=psi)
  pp = sg.Pointset(pp)
  vv1 = sg.eval_layer(G, pp)
  vv2 = ly.fundsol(abs(pp.x - z0), 0)

  h = plt.figure()
  plot.plot(h, x, y, vv2 + vv1,'cf')
  plt.show(block=True)
