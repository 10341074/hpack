#!../../venv/bin/python
import matplotlib.pyplot as plt
import numpy as np

from scipy.linalg import solve
import scipy.linalg as linalg
import numpy.linalg

from __types__ import *
import layerpot as ly
import segment as sg
import plot

verbose=0

# (K' + 0.5 * c(h) * I) psi = -phi_nu

def directA(l, c): # to delete
  Kp = l.exp(0, l.b[0], [], [])
  A = Kp
  A[np.diag_indices(len(l.b[0].x))]=A[np.diag_indices(len(l.b[0].x))] + 0.5 * c
  return A

def directKpc(l, c):
  # Kp = l.exp(s=l.b[0])
  Kp = l.eval_self(exp=ly.layerpotSD)
  Kpc = Kp
  Kpc[np.diag_indices(l.n)] = Kpc[np.diag_indices(l.n)] + 0.5 * c
  if verbose:
    print('directKpc condition number= ', numpy.linalg.cond(np.array(Kpc, float)))
  return Kpc

def directrhs(l=[], b=[], z0=[]):
  if z0 == []:
    print('Error: z0 in directrhs')
  if l != []:
    d = np.array([z for bk in l.b for z in bk.x]) - z0
    nx = np.array([nxk for bk in l.b for nxk in bk.nx])
  elif b != []:
    d = b.x - z0
    nx = b.nx
  r = abs(d)
  cosphi = ly.scalar(nx, d) / r
  rhs = - - ly.fundsol_deriv(r, cosphi, 0)
  return rhs
  
def directpb(l, c, z0, A_f=((), ()), rhs=[]):
  A, solutor = A_f
  if A == ():
    print('Warning: calculated A, not passed')
    A = directKpc(l, c)
    solutor = linalg.solve
  if rhs == []:
    rhs = directrhs(l=l, z0=z0)
  psi = solutor(A, rhs)
  return psi

def plotdpb(l, z0, x1_x2_xn, y1_y2_yn=((), (), ()), psi=(), t='im'):
  x, y, pp = plot.meshgrid(x1_x2_xn, y1_y2_yn)
  lPsi = sg.Layer(b=l.b, exp=ly.layerpotS, dns=l.dns)
  if psi != ():
    lPsi.dns = psi
  pp = sg.Pointset(pp)
  #uPsi = sg.eval_layer(lPsi, pp)
  uPsi = ly.layerpotS(s=l.b[0], t=pp).dot(lPsi.dns)
  uPhi = ly.fundsol(abs(pp.x - z0), 0)
  plot.plot(x, y, uPsi + uPhi, t=t, show=0)
  l.plot(p=True)
  # so.plot(p=True)
  # sb.plot(p=True)
  # print('z0', z0 )
  plt.show(block=False)
  return

if __name__ == "__main__":
  h = 500
  n = 50
  c = 1. * (h + 1)/(h-1)

  # (K' + 0.5 * c(h) * I) psi = -phi_nu
  sd = sg.Segment(100, Z=sg.dZ, Zp=sg.dZp, Zpp=sg.dZpp, args=[], quad='gp')
  s = sg.Segment(100, f = sg.circle, inargs = (-1, 0.5), quad='ps')
  # s = sg.Segment(n, Z=sg.kZ, Zp=sg.kZp, Zpp=sg.kZpp, args=[], periodic=True)
  b = sg.Boundary([sd])
  b2 = sg.Boundary([s])
  l = sg.Layer([b, b2], ly.layerpotSD)
  l1 = sg.Layer([b], ly.layerpotSD)
  l2 = sg.Layer([b2], ly.layerpotSD)

  z0 = 5.1
  d = s.x - z0
  r = abs(d)
  cosphi = np.real(np.conj(s.nx) * d) / r
  # rhs = - ly.fundsol_deriv(r, cosphi, 0)
  Kpc = directKpc(l=l, c=c)
  rhs = directrhs(l=l, z0=z0)
  Kpc1 = directKpc(l=l1, c=c)
  rhs1 = directrhs(l=l1, z0=z0)
  Kpc2 = directKpc(l=l2, c=c)
  rhs2 = directrhs(l=l2, z0=z0)

  psi = solve(Kpc, rhs)
  p1 = solve(Kpc1, rhs1)
  p2 = solve(Kpc2, rhs2)
  print(p1.shape)
  print(p2.shape)
  
  psi2 = np.concatenate((p1, p2))
  x, y, pp = plot.meshgrid((-1.3, 1.3, 1200))
  x, y, pp = plot.meshgrid((-2, 4, 80),(-2, 2, 60))
  G = sg.Layer(b=[b, b2], exp=ly.layerpotS, dns=psi)
  pp = sg.Pointset(pp)
  vv1 = sg.eval_layer(G, pp)
  vv2 = ly.fundsol(abs(pp.x - z0), 0)

  plot.plot(x, y, vv2 + vv1, 'srf')
  # s.plot(p=True)
  plt.show(block=True)
  
def mapNtoD(l, g):
  Kp = l.eval_self(exp=ly.layerpotSD)
  Kp2 = ly.layerpotSD(s=l.b[0])
  # print('check ', max([max(abs(r)) for r in Kp - Kp2]))
  
  Kpa = Kp
  Kpa[np.diag_indices(l.n)] = Kpa[np.diag_indices(l.n)] - 0.5
  if verbose:
    print('mapNtoD condition number= ', numpy.linalg.cond(np.array(Kpa, float)))
  phi = linalg.solve(Kpa, g)
  sl = sg.Layer(b=l.b, exp=ly.layerpotS, dns=phi)
  S = sl.eval_self()
  trace = S.dot(sl.dns)
  return trace
  
