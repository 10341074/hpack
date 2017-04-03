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

verbose = 1

def computeRtikh(R, a):
  return a * np.eye(R.shape[1]) + R.T.dot(R)

def computeRHS(R, U, U_nu, z0, so, theta=0, RHS=[]):
  if RHS != []:
    return RHS
  F = np.cos(theta) * ly.phi_x(z0, so.x) + np.sin(theta) * ly.phi_y(z0, so.x)
  F_nu = np.cos(theta) * ly.phi_x_n(z0, so.x, so.nx) + np.sin(theta) * ly.phi_y_n(z0, so.x, so.nx)
  return U.dot(F_nu) - U_nu.dot(F)

def computeRHStikh(R, U, U_nu, z0, so, theta=0, RHS=[]):
  if RHS == []:
    RHS = computeRHS(R, U, U_nu, z0, so, theta)
  return R.T.dot(RHS)

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
    V2_nu = ly.layerpotDD(s=sb, t=so)
    R2 = U.dot(V2_nu) - U_nu.dot(V2)
    R = R2
  if testset == 3:
    R = np.concatenate((R1.T, R2.T)).T
  if verbose:
    print('R.shape = ', R.shape)
  return (R, U, U_nu)

def gap_init(R, a, reg, regmet, solver):
  if reg or reg ==1:
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
  ninv = np.empty(len(pp.x), float)
  #res = np.empty((len(pp.x), 1), float)
  #nsolgap= np.empty((len(pp.x), 1), float)
  res=()
  nsolgap=()
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
    # if abs(z0) < 1e-10:
    #   print('k= ',computek(zeta, z0, ld2, lb))
    #nn[k] = computek(zeta, z0, ld, lb)
    #normf(k) = real(z0) * imag(z0) + real(z0);
    #normf[k] = max(abs(solve(Rreg, RHSreg) - zeta))
    # print('residual= ',linalg.norm(Rreg.dot(zeta)-RHSreg))
    # print('residual2= ',linalg.norm(R.dot(zeta)-RHS))
    #sum(F_nu .* so.w.'); % check null
  return (ninv, res, nsolgap)

def meshgrid(x1_x2_xn, y1_y2_yn=((), (), ())):
  x, y, pp = plot.meshgrid(x1_x2_xn, y1_y2_yn)
  pp = sg.Pointset(pp)
  return (x, y, pp)
