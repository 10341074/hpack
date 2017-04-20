# from domain import *
import matplotlib.pyplot as plt
import numpy as np
# from scipy.linalg import solve
# from scipy.linalg import lstsq
from numpy.linalg import norm
from numpy.linalg import cond
import scipy.linalg as linalg
import numpy.linalg
import time

from __types__ import *
import shapes as sh
import layerpot as ly
import segment as sg
import plot
import directproblem as dpb

import linfunc as linf
import lintype as lint

verbose = 1

def computeRtikh(R, a):
  return a * np.eye(R.shape[1]) + R.T.dot(R)

def computeRHS(R, U, U_nu, z0, so, theta=0, RHS=[]):
  if RHS != []:
    return RHS
  F = np.cos(theta) * ly.phi_x(z0, so.x) + np.sin(theta) * ly.phi_y(z0, so.x)
  F_nu = np.cos(theta) * ly.phi_x_n(z0, so.x, so.nx) + np.sin(theta) * ly.phi_y_n(z0, so.x, so.nx)
  return U.dot(F_nu) - U_nu.dot(F)

def computeRHStikhold(R, U=(), U_nu=(), z0=(), so=(), theta=0, RHS=[]):
  if RHS == []:
    RHS = computeRHS(R, U, U_nu, z0, so, theta)
  return R.T.dot(RHS)
def computeRHStikh(isolver, RHS=()):
  if RHS == ():
    print('Error, RHStikh cannot compute RHS')
  A = isolver.ARHS
  return A.T.dot(RHS)

def computeAtikh(A, alpha):
  return alpha * np.eye(A.shape[1]) + A.T.dot(A)

def computeRlstsq(R):
  return R.T.dot(R)

def computeRHSlstsq(R, U, U_nu, z0, so, RHS=[]):
  if RHS == []:
    RHS = computeRHS(R, U, U_nu, z0, so)
  return R.T.dot(RHS)

def computek(zeta, z0, ld, lb):
  F = ly.phi(z0, ld.x)
  F_nu = ly.phi_n(z0, ld.x, ld.nx)
  APsi = dpb.mapNtoD(ld, F_nu)
  Psi = F
  vz = ly.layerpotS(s=lb, t=ld)
  vz = vz.dot(zeta)
  i_v = (vz**2).dot(ld.w)
  i_A = (APsi**2).dot(ld.w)
  i_P = (Psi**2).dot(ld.w)
  i_AP = (APsi * Psi).dot(ld.w)
  # (AP - v + (AP^2 - 2*v*AP + A*v + P*v - A*P)^(1/2))/(P - v)
  # k = (i_AP - i_v + (i_AP**2 - 2*i_v*i_AP + i_A*i_v + i_P*i_v - i_A*i_P)**(0.5))/(i_P - i_v)
  nn = (1 - h)*vz - APsi + h * Psi
  nn = (nn**2).dot(ld.w)
  return nn

def computeallpsi(ld, sb, c):
  nsd = ld.n
  nsb = sb.n
  allpsi = np.empty((nsd, nsb), float)
  A = dpb.directKpc(ld, c)
  for k in range(nsb):
    z0 = sb.x[k]
    rhs = dpb.directrhs(l=ld, z0=z0)
    allpsi[:, k] = dpb.directpb(l=ld, c=c, z0=z0, A_f = (A, linalg.solve), rhs=rhs)
  return allpsi

def computeR(allpsi, ld, so, sb, sv=(), testset=1):
  if sv == ():
    sv = sb
  sd = ld.b[0]
  nso, nsb = so.n, sb.n

  U_psi = np.empty((nso, nsb), float)
  U_psi_nu = np.empty((nso, nsb), float)
  # kerS = sg.eval_layer(ld, so, exp = ly.layerpotS)
  # kerSD = sg.eval_layer(ld, so, exp = ly.layerpotSD)
  kerS = ly.layerpotS(0, sd, so)
  kerSD = ly.layerpotSD(0, sd, so)

  for k in range(nsb):
    U_psi[:,k] = kerS.dot(allpsi[:,k])
    U_psi_nu[:,k] = kerSD.dot(allpsi[:,k])

  # U_psi = np.diag(so.w) * U_psi
  # U_psi_nu = np.diag(so.w) * U_psi_nu
  U_psi = U_psi.T.dot(np.diag(so.w))
  U_psi_nu = U_psi_nu.T.dot(np.diag(so.w))

  U = ly.layerpotS(s=so, t=sb)
  U_nu = ly.layerpotD(s=so, t=sb)

  # U = U + U_psi.';
  # U_nu = U_nu + U_psi_nu.';
  
  U = U + U_psi
  U_nu = U_nu + U_psi_nu

  #lb = sg.Layer([sd], ly.layerpotS)
  testset = 1
  if testset == 1 or testset == 3:
    V1 = ly.layerpotS(s=sv, t=so)
    V1_nu = ly.layerpotSD(s=sv, t=so)
    R1 = U.dot(V1_nu) - U_nu.dot(V1)
    R = R1
  if testset == 2 or testset == 3:
    V2 = ly.layerpotD(s=sv, t=so)
    V2_nu = ly.layerpotDD(s=sv, t=so)
    R2 = U.dot(V2_nu) - U_nu.dot(V2)
    R = R2
  if testset == 3:
    R = np.concatenate((R1.T, R2.T)).T
  if verbose:
    print('R.shape = ', R.shape)
  return (R, U, U_nu)

def gap_init(R, a, reg, regmet, solver):
  if reg or reg == 1:
    if regmet == 'tikh':
      Rreg = computeRtikh(R, a)
      _A_b = (Rreg, computeRHStikh)
  else:
    _A_b = (R, computeRHS)

  if solver == 's':
    # _gap_s = (linalg.solve, _A_b[0], _A_b[1])
    _gap_s = (_solve, _A_b[0], _A_b[1])
    _gap = _gap_s
  elif solver == 'lu':
    (lu, piv) = linalg.lu_factor(_A_b[0])
    # _gap_lu = (linalg.lu_solve, (lu, piv), _A_b[1])
    _gap_lu = (_lu, (lu, piv), _A_b[1])
    _gap = _gap_lu
  elif solver == 'lstsq':
    # _gap_lstsq = (linalg.lstsq, _A_b[0], _A_b[1])
    _gap_lstsq = (_lstsq, _A_b[0], _A_b[1])
    _gap = _gap_lstsq # to change
  else:
    print('Error: not valid solver')
  return _gap

def _solve(_m, _rhs):
  return linalg.solve(_m, _rhs)
def _lu(_m, _rhs):
  return linalg.lu_solve(_m, _rhs)
def _lstsq(_m, _rhs):
  return linalg.lstsq(_m, _rhs)[0]


def computeallsolsgap(_gap, pp, R, U, U_nu, so, theta=0):
  res = ()
  nsolgap = ()
  ninv = np.empty(len(pp.x), float)
  res = np.empty((len(pp.x), 1), float)
  # nsolgap= np.empty((len(pp.x), 1), float)
  # nn = np.empty(len(pp.x), float)
  # (lu, piv) = linalg.lu_factor(computeRtikh(R, 1e-14))
  for k in range(len(pp.x)):
    z0 = pp.x[k]
    RHS = computeRHS(R, U, U_nu, z0, so, theta)
    # RHSreg = computeRHStikh(R, U, U_nu, z0, so, theta, RHS=RHS)
    _gap_rhs =  _gap[2](R, U, U_nu, z0, so, theta, RHS=RHS)
    zeta = _gap[0](_gap[1], _gap_rhs)
    # zeta = _gap[0](_gap[1], _gap[2](R, U, U_nu, z, so, theta))[0]
    # zeta = linalg.lu_solve((lu, piv), RHSreg)
    # res[k]= norm(R.dot(zeta) - RHS)
    # nsolgap[k]= norm(zeta)
    time.sleep(0.001)  
    ninv[k] = norm(RHS) / norm(zeta)
    res[k] = linalg.norm(R.dot(zeta) - RHS)
    # print('residual= ',linalg.norm(Rreg.dot(zeta)-RHSreg))
    # print('residual2= ',linalg.norm(R.dot(zeta)-RHS))
  return (ninv, res, nsolgap)

def meshgrid(x1_x2_xn, y1_y2_yn=((), (), ())):
  x, y, pp = plot.meshgrid(x1_x2_xn, y1_y2_yn)
  pp = sg.Pointset(pp)
  return (x, y, pp)

def computeallpsiL0(so, T):
  # L = np.empty((len(T), so.n))
  allpsi0 = dpb.mapNtoD0(so, T, so.s0) 
  Lo = ly.layerpotS(s=so)
  L = Lo.dot(allpsi0)
  means = np.ones(so.n).dot(np.diagflat(so.w).dot(L)) / sum(so.w)
  L = L - np.array([means for k in range(so.n)])
  return L
def computeallpsiL(so, ld, T, c):
  # L = np.empty((len(T), so.n))
  allpsi = dpb.mapNtoD(so, ld, T, c, so.s0)
  Lo = ly.layerpotS(s=so)
  Ld = ly.layerpotS(s=ld, t=so)
  L = Lo.dot(allpsi[0:so.n]) + Ld.dot(allpsi[so.n::])
  means = sum(np.diagflat(so.w).dot(L)) / sum(so.w)
  L = L - np.array([means for k in range(so.n)])
  return T.T.dot(L)

def computeL0(so, T):
  allpsi0 = dpb.mapNtoD0(so, T, so.s0) 
  Lo = ly.layerpotS(s=so)
  L0 = Lo.dot(allpsi0)
  means = np.ones(so.n).dot(np.diagflat(so.w).dot(L0)) / sum(so.w)
  L0 = L0 - np.array([means for k in range(so.n)])
  return L0
def computeL0B(so):
  L0 = computeL0(so, so.B[:, 1:])
  return so.B[:, 1:].T.dot(np.diagflat(so.w).dot(L0))
def computeL(so, ld, T, c):
  allpsi = dpb.mapNtoD(so, ld, T, c, so.s0)
  Lo = ly.layerpotS(s=so)
  Ld = ly.layerpotS(s=ld, t=so)
  L = Lo.dot(allpsi[0:so.n]) + Ld.dot(allpsi[so.n::])
  # means = sum(np.diagflat(so.w).dot(L)) / sum(so.w) # correct? strange sum by rows
  means = np.ones(so.n).dot(np.diagflat(so.w).dot(L)) / sum(so.w)
  L = L - np.array([means for k in range(so.n)])
  return L
def computeLB(so, ld, c):
  L = computeL(so, ld, so.B[:, 1:], c)
  return so.B[:, 1:].T.dot(np.diagflat(so.w).dot(L))

def computeLL0(ld, so, sb, c, testset=0):
  # Kp = ly.layerpotSD(s=so)
  # (nu, s0) = linf.eigmaxpowerw(A=Kp, s=so)
  # if testset == 0:
  #   S = linf.gramschmidtw(so, s0=s0)
  #   T = S[:, 1::]
  # if testset == 1 or testset == 3:
  #   V1_nu = ly.layerpotSD(s=sb, t=so)
  #   T = V1_nu
  # if testset == 2 or testset == 3:
  #   V2_nu = ly.layerpotDD(s=sb, t=so)
  #   T = V2_nu
  # if testset == 3:
  #   T = np.concatenate((V1_nu.T, V2_nu.T)).T

  L0 = computeL0(so, so.B[:, 1:])
  L = computeL(so, ld, so.B[:, 1:], c)

  LL0 = L - L0
  return LL0
def computeLL0B(ld, so, sb, c, testset=0):
  L0 = computeL0B(so)
  L = computeLB(so, ld, c)

  LL0 = L - L0
  return LL0

def computeRHSNtoD(LL0, z0, so, theta=0, RHS=()):
  if RHS != ():
    return RHS
  S = linf.gramschmidtw(so)
  RHS = S[:, 1::].T.dot(ly.phi_theta(z0, so.x, theta))
  return RHS

def computeRHSNtoDtikh(LL0, z0, so, theta=0, RHS=()):
  if RHS == ():
    RHS = computeRHSNtoD(LL0, z0, so, theta)
  return LL0.T.dot(RHS)

def NtoD_init(LL0, a, reg, regmet, solver):
  if reg or reg == 1:
    if regmet == 'tikh':
      Rreg = computeRtikh(LL0, a)
      _A_b = (Rreg, computeRHSNtoDtikh)
  else:
    _A_b = (LL0, computeRHSNtoD)

  if solver == 's':
    # _gap_s = (linalg.solve, _A_b[0], _A_b[1])
    _gap_s = (_solve, _A_b[0], _A_b[1])
    _gap = _gap_s
  elif solver == 'lu':
    (lu, piv) = linalg.lu_factor(_A_b[0])
    # _gap_lu = (linalg.lu_solve, (lu, piv), _A_b[1])
    _gap_lu = (_lu, (lu, piv), _A_b[1])
    _gap = _gap_lu
  elif solver == 'lstsq':
    # _gap_lstsq = (linalg.lstsq, _A_b[0], _A_b[1])
    _gap_lstsq = (_lstsq, _A_b[0], _A_b[1])
    _gap = _gap_lstsq # to change
  else:
    print('Error: not valid solver')
  return _gap

def solver_init(A, alpha, reg, regmet, solver, RHS_fcom, RHS_args):
  s = lint.Solver(A=A, RHS_fcom=RHS_fcom, RHS_args=RHS_args)
  if reg or reg == 1:
    if regmet == 'tikh':
      s.A = computeAtikh(A, alpha)
      s.ARHS = A
      s.RHS_f = computeRHStikh
    else:
      print('regularization without specified method')
  else:
    s.A = A
    s.ARHS = A
    s.RHS_f = s.RHS_fcom

  if solver == 's':
    s.solver_a = s.A
    s.solver_f = _solve
  elif solver == 'lu':
    (lu, piv) = linalg.lu_factor(s.A)
    s.solver_a = (lu, piv)
    s.solver_f = _lu
  elif solver == 'lstsq':
    s.solver_a = s.A
    s.solver_f = _lstsq # to change
  else:
    print('Error: not valid solver')
  return s

def NtoD_computeneumb(args):
  s, z0, theta = args['s'], args['z0'], args['theta']
  xx = ly.phi_xx(z0, s.x)
  xy = ly.phi_xy(z0, s.x)
  yy = ly.phi_yy(z0, s.x)
  # hess = [[xx, xy] , [xy, yy]]
  a = np.cos(theta) + 1j * np.sin(theta)
  rhs = ly.scalar(a, ly.phi_p(z0, s.x))
  neum = s.nx.real * xx * a.real + s.nx.real * xy * a.imag + \
         s.nx.imag * xy * a.real + s.nx.imag * yy * a.imag
  neumb = neum.dot(np.diagflat(s.w).dot(s.B))
  return neumb
def NtoD_computeRHS(args, rhs=()):
  L0, L0B, s, z0, theta = args['L0'], args['L0B'], args['s'], args['z0'], args['theta']
  if rhs != (): # remove check at every call
    return rhs
  a = np.cos(theta) + 1j * np.sin(theta)
  neumb = NtoD_computeneumb(args)
  dirh = L0.dot(neumb[1:])
  rhs = ly.scalar(a, ly.phi_p(z0, s.x))
  m = sum(rhs * s.w) / sum(s.w)
  rhs = rhs - m
  return rhs - dirh

def NtoD_computeRHSB(args, rhs=()):
  L0, L0B, s, z0, theta = args['L0'], args['L0B'], args['s'], args['z0'], args['theta']
  if rhs != (): # remove check at every call
    return rhs
  a = np.cos(theta) + 1j * np.sin(theta)
  # rhs = NtoD_computeRHS(args)
  neumb = NtoD_computeneumb(args)
  dirh = L0B.dot(neumb[1:])
  rhs = ly.scalar(a, ly.phi_p(z0, s.x))
  m = sum(rhs * s.w) / sum(s.w)
  rhs = rhs - m
  rhs = s.B[:, 1:].T.dot(np.diagflat(s.w).dot(rhs))
  return rhs - dirh
def rhs(z0, so, theta=0, s0=()):
  # S = linf.gramschmidtw(so)
  # r1 = S[:, 1::].T.dot(ly.phi_theta(z0, so.x, theta))
  r1  = ly.phi_theta(z0, so.x, theta)
  xx = ly.phi_xx(z0, so.x)
  xy = ly.phi_xy(z0, so.x)
  yy = ly.phi_yy(z0, so.x)
  g = so.nx.real * xx * np.cos(theta) + so.nx.real * xy * np.sin(theta) + so.nx.imag * xy * np.cos(theta) + so.nx.imag * yy * np.sin(theta)
  # psi = dpb.mapNtoD0(so, -g, so.s0)
  # f = r1 + ly.layerpotS(s=so).dot(psi)
  # f = f - sum(f*so.w) / sum(so.w) * np.ones(so.n)
  return g

def computeallsolsNtoD(_NtoD, pp, LL0, so, theta=0):
  res = ()
  nsolgap = ()
  ninv = np.empty(len(pp.x), float)
  S = linf.gramschmidtw(so)
  s0 = linf.eigmaxpowerw(ly.layerpotS(s=so), so)[1]
  # res = np.empty((len(pp.x), 1), float)
  # nsolgap= np.empty((len(pp.x), 1), float)
  # nn = np.empty(len(pp.x), float)
  # (lu, piv) = linalg.lu_factor(computeRtikh(R, 1e-14))
  for k in range(len(pp.x)):
    z0 = pp.x[k]
    # RHS = S[:, 1::].T.dot(ly.phi_theta(z0, so.x, theta))
    # RHS = S[:, 1::].T.dot(rhs(z0, so, theta, s0))
    RHS =rhs(z0, so, theta, s0)
    # RHS = ly.phi_theta(z0, so.x, theta)
    #RHSreg = LL0.T.dot(RHS)
    _NtoD_rhs =  _NtoD[2](LL0, z0, so, theta, RHS=RHS)
    zeta = _NtoD[0](_NtoD[1], _NtoD_rhs)
    # zeta = _gap[0](_gap[1], _gap[2](R, U, U_nu, z, so, theta))[0]
    # zeta = linalg.lu_solve((lu, piv), RHSreg)
    # res[k]= norm(R.dot(zeta) - RHS)
    # nsolgap[k]= norm(zeta)
    time.sleep(0.01)  
    ninv[k] = norm(RHS) / norm(zeta)
    # print('residual= ',linalg.norm(Rreg.dot(zeta)-RHSreg))
    # print('residual2= ',linalg.norm(R.dot(zeta)-RHS))
  return (ninv, res, nsolgap)

def iallsols(isolver, pointstest, LL0, so, theta=0):
  res = ()
  nsolgap = ()
  ninv = np.empty(len(pointstest.x), float)
  # S = linf.gramschmidtw(so)
  # res = np.empty((len(pp.x), 1), float)
  # nsolgap= np.empty((len(pp.x), 1), float)
  # nn = np.empty(len(pp.x), float)
  # (lu, piv) = linalg.lu_factor(computeRtikh(R, 1e-14))
  for k in range(len(pointstest.x)):
    isolver.RHS_args['z0'] = pointstest.x[k]
    # RHS = S[:, 1::].T.dot(ly.phi_theta(z0, so.x, theta))
    # RHS = ly.phi_theta(z0, so.x, theta)
    #RHSreg = LL0.T.dot(RHS)
    RHS = isolver.RHS_fcom(isolver.RHS_args)
    RHSr =  isolver.RHS_f(isolver, RHS=RHS)
    zeta = isolver.solver_f(isolver.solver_a, RHSr)
    # zeta = _gap[0](_gap[1], _gap[2](R, U, U_nu, z, so, theta))[0]
    # zeta = linalg.lu_solve((lu, piv), RHSreg)
    # res[k]= norm(R.dot(zeta) - RHS)
    # nsolgap[k]= norm(zeta)
    time.sleep(0.001)  
    ninv[k] = linalg.norm(RHS) / linalg.norm(zeta)
    # print('residual= ',linalg.norm(Rreg.dot(zeta)-RHSreg))
    # print('residual2= ',linalg.norm(R.dot(zeta)-RHS))
  return (ninv, res, nsolgap)

