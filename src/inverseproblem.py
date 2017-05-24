# from domain import *
import matplotlib.pyplot as plt
import numpy as np
from numpy.linalg import norm
from numpy.linalg import cond
import scipy.linalg as linalg
import numpy.linalg
import scipy.stats
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

def gap_computeRHS(args, RHS=()):
  R, U, U_nu, z0, so, theta = args['R'], args['U'], args['U_nu'], args['z0'], args['so'], args['theta']
  if RHS != ():
    return RHS
  F = np.cos(theta) * ly.phi_x(z0, so.x) + np.sin(theta) * ly.phi_y(z0, so.x)
  F_nu = np.cos(theta) * ly.phi_x_n(z0, so.x, so.nx) + np.sin(theta) * ly.phi_y_n(z0, so.x, so.nx)
  return U.dot(F_nu) - U_nu.dot(F)
def gap_computeRHSB(args, RHS=()):
  RHS = gap_computeRHS(args, RHS)
  return s.B[:, 1:].T.dot(np.diagflat(s.w).dot(RHS))

def computeRHStikh(isolver, RHS=()):
  if RHS == ():
    print('Error, RHStikh cannot compute RHS')
  A = isolver.A
  return A.T.dot(RHS)

def computeAtikh(A, alpha):
  return alpha * np.eye(A.shape[1]) + A.T.dot(A)

# def computek(zeta, z0, ld, lb):
#   F = ly.phi(z0, ld.x)
#   F_nu = ly.phi_n(z0, ld.x, ld.nx)
#   APsi = dpb.mapNtoD(ld, F_nu)
#   Psi = F
#   vz = ly.layerpotS(s=lb, t=ld)
#   vz = vz.dot(zeta)
#   i_v = (vz**2).dot(ld.w)
#   i_A = (APsi**2).dot(ld.w)
#   i_P = (Psi**2).dot(ld.w)
#   i_AP = (APsi * Psi).dot(ld.w)
#   # (AP - v + (AP^2 - 2*v*AP + A*v + P*v - A*P)^(1/2))/(P - v)
#   # k = (i_AP - i_v + (i_AP**2 - 2*i_v*i_AP + i_A*i_v + i_P*i_v - i_A*i_P)**(0.5))/(i_P - i_v)
#   nn = (1 - h)*vz - APsi + h * Psi
#   nn = (nn**2).dot(ld.w)
#   return nn

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

def computeU(allpsi, ld, so, sb):
  nso, nsb = so.n, sb.n
  U_psi = np.empty((nso, nsb), float)
  U_psi_nu = np.empty((nso, nsb), float)
  # kerS = sg.eval_layer(ld, so, exp = ly.layerpotS)
  # kerSD = sg.eval_layer(ld, so, exp = ly.layerpotSD)
  kerS = ly.layerpotS(0, ld, so)
  kerSD = ly.layerpotSD(0, ld, so)

  # for k in range(nsb):
  #   U_psi[:,k] = kerS.dot(allpsi[:,k])
  #   U_psi_nu[:,k] = kerSD.dot(allpsi[:,k])
  U_psi = kerS.dot(allpsi)
  U_psi_nu = kerSD.dot(allpsi)

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
  return (U, U_nu)

def computeR(allpsi, ld, so, sb, sv=(), testset=1):
  if sv == ():
    sv = sb
  nso, nsb = so.n, sb.n

  U, U_nu = computeU(allpsi, ld, so, sb)
  
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

def _solve(_m, _rhs):
  return linalg.solve(_m, _rhs)
def _lu(_m, _rhs):
  return linalg.lu_solve(_m, _rhs)
def _lstsq(_m, _rhs):
  return linalg.lstsq(_m, _rhs)[0]

def meshgrid(x1_x2_xn, y1_y2_yn=((), (), ())):
  x, y, pp = plot.meshgrid(x1_x2_xn, y1_y2_yn)
  pp = sg.Pointset(pp)
  return (x, y, pp)
#####################
def computeL0(so, T):
  if T == ():
    T = so.BX
  print('computing L0')
  allpsi0 = dpb.mapNtoD0(so, T, so.s0) 
  Lo = ly.layerpotS(s=so)
  L0 = Lo.dot(allpsi0)
  means = np.ones(so.n).dot(np.diagflat(so.w).dot(L0)) / sum(so.w)
  L0 = L0 - np.array([means for k in range(so.n)])
  return L0
def computeL0B(so, T, L0=()):
  if L0 == ():
    L0 = computeL0(so, T)
    print('computing L0B')
  else:
    print('computing L0B from L0')
  # return so.BY.T.dot(np.diagflat(so.w).dot(L0))
  return so.BYinv.dot(L0)
def computeL(ld, so, T, c):
  if T == ():
    T = so.BX
  print('computing L')
  allpsi = dpb.mapNtoD(so, ld, T, c, so.s0)
  Lo = ly.layerpotS(s=so)
  Ld = ly.layerpotS(s=ld, t=so)
  L = Lo.dot(allpsi[0:so.n]) + Ld.dot(allpsi[so.n::])
  # means = sum(np.diagflat(so.w).dot(L)) / sum(so.w) # correct? strange sum by rows
  means = np.ones(so.n).dot(np.diagflat(so.w).dot(L)) / sum(so.w)
  L = L - np.array([means for k in range(so.n)])
  return L
def computeLB(ld, so, T, c, L=()):
  if L == ():
    L = computeL(ld, so, T, c)
    print('computing LB')
  else:
    print('computing LB from L')
  # return so.BY.T.dot(np.diagflat(so.w).dot(L))
  return so.BYinv.dot(L)
def computeLL0(ld, so, T, c, L0=(), L=()):
  if T == ():
    T = so.BX
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
  if L0 == ():
    L0 = computeL0(so, T)
    print('computing L0/LL0')
  else:
    print('taking L0/LL0 from L0')

  if L == ():
    L = computeL(ld, so, T, c)
    print('computing L/LL0')
  else:
    print('taking L/LL0 from L')

  LL0 = L - L0
  return LL0
def computeLL0B(ld, so, T, c, L0B=(), LB=()):
  if L0B == ():
    L0B = computeL0B(so, T)
    print('computing L0B/LL0B')
  else:
    print('taking L0B/LL0B from L0B')
  if LB == ():
    LB = computeLB(ld, so, T, c)
    print('computing LB/LL0B')
  else:
    print('taking LB/LL0B from LB')

  LL0 = LB - L0B
  return LL0
# def computeLL0B(ld, so, T, c, LL0=():
#   if LL0 == ():
#     LL0 = computeL(ld, so, T, c)
#     print('computing LB')
#   else:
#     print('computing LB from L')
#   return 
def computeLLdiff(ld, so, T, c):
  if T == ():
    T = so.BX
  allpsi0 = dpb.mapNtoD0(so, T, so.s0) 
  Lo = ly.layerpotSD(s=so, t =ld)
  rhsdiff = Lo.dot(allpsi0)
  allpsi = dpb.mapNtoDdiff(so, ld, rhsdiff, c, so.s0)
  Lo = ly.layerpotS(s=so)
  Ld = ly.layerpotS(s=ld, t=so)
  L = Lo.dot(allpsi[0:so.n]) + Ld.dot(allpsi[so.n::])
  # means = sum(np.diagflat(so.w).dot(L)) / sum(so.w) # correct? strange sum by rows
  means = np.ones(so.n).dot(np.diagflat(so.w).dot(L)) / sum(so.w)
  L = L - np.array([means for k in range(so.n)])
  return L
def computeLLBdiff(ld, so, T, c):
  if T == ():
    T = so.B[:, 1:]
  L = computeLLdiff(ld, so, T, c)
  return so.BY.T.dot(np.diagflat(so.w).dot(L))
  
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

def solver_init(A, alpha, reg, regmet, solver, RHS_fcom, RHS_args, BX=(), BY=()):
  s = lint.Solver(A=A, RHS_fcom=RHS_fcom, RHS_args=RHS_args, BX=BX, BY=BY)
  s.alpha = alpha
  s.delta = 1e-6
  if reg or reg == 1:
    if regmet == 'tikh':
      s.Ar = computeAtikh(A, alpha)
      s.RHS_f = computeRHStikh
    else:
      print('regularization without specified method')
  else:
    s.RHS_f = s.RHS_fcom

  if solver == 's':
    s.solver_a = s.Ar
    s.solver_f = _solve
  elif solver == 'lu':
    (lu, piv) = linalg.lu_factor(s.Ar)
    s.solver_a = (lu, piv)
    s.solver_f = _lu
  elif solver == 'lstsq':
    s.solver_a = s.Ar
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
  # rhs = ly.scalar(a, ly.phi_p(z0, s.x))
  neum = s.nx.real * xx * a.real + s.nx.real * xy * a.imag + \
         s.nx.imag * xy * a.real + s.nx.imag * yy * a.imag
  # modified !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  # neumb = neum.dot(np.diagflat(s.w).dot(s.B))
  return neum

def NtoD_computeRHS(args, rhs=()):
  L0, L0B, s, z0, theta = args['L0'], args['L0B'], args['s'], args['z0'], args['theta']
  if rhs != (): # remove check at every call
    return rhs
  a = np.cos(theta) + 1j * np.sin(theta)
  neumb = NtoD_computeneumb(args)
  # print('mean', sum(neumb*s.w))
  # need decomposition BX of neumb !!!!!!!!!!!!!!!!!!!
  dirh = L0.dot(neumb)
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
  dirh = L0B.dot(neumb)
  rhs = ly.scalar(a, ly.phi_p(z0, s.x))
  m = sum(rhs * s.w) / sum(s.w)
  rhs = rhs - m
  rhs = s.BY.T.dot(np.diagflat(s.w).dot(rhs))
  return rhs - dirh

def func_disc_pow1_n(res, sol, w, delta):
  # \|Kx^{\alpha,\delta} - y^\delta\| - \delta \|x^{\alpha,\delta}\|
  return np.sqrt(sum(res**2 * w)) - delta * np.sqrt(sum(sol**2 * w))
def func_disc_p_pow1_n(res, Kdx, sol, sol_p, w, delta):
  # \frac{1}{\|Kx^{\alpha,\delta} - y^\delta\|} \langle  Kx^{\delta, \alpha} - y^\delta, K\frac{d}{d\alpha}x^{\alpha, \delta}\rangle -\delta\frac{1}{\|x^{\alpha,\delta}\|}\langle x^{\alpha, \delta}, \frac{d}{d\alpha}x^{\alpha, \delta}\rangle
  return sum(res * Kdx * w) / np.sqrt(sum(res**2 * w)) - delta * sum(sol * sol_p * w) / np.sqrt(sum(sol**2 *w))


def func_alpha_newton(alpha, disc, disc_p):
  return alpha - disc / disc_p

def func_alpha_bis(alpha, alpha_l, alpha_r, disc):
  if disc > 0:
    alpha_r = alpha
  else:
    alpha_l = alpha
  alpha = (alpha_l + alpha_r) * 0.5
  # print('from bisec: alpha= ',alpha)
  return (alpha, alpha_l, alpha_r)

def iallsols(isolver, pointstest, so):
  w = so.w
  isolver.save_zeta = np.empty((len(pointstest.x), isolver.A.shape[1]), float)
  isolver.save_sol = np.empty((len(pointstest.x), so.n), float)
  isolver.save_ratio = np.empty((len(pointstest.x),1), float)

  # test one point
  isolver.RHS_args['z0'] = pointstest.x[0]
  RHS = isolver.RHS_fcom(isolver.RHS_args)

  isolver.save_rhs = np.empty((len(pointstest.x), len(RHS)), float)
  for k in range(len(pointstest.x)):
    isolver.RHS_args['z0'] = pointstest.x[k]
    RHS = isolver.RHS_fcom(isolver.RHS_args)
    RHSr = isolver.RHS_f(isolver, RHS=RHS)
    zeta = isolver.solver_f(isolver.solver_a, RHSr)
    isolver.save_zeta[k, :] = zeta
    isolver.save_rhs[k] = RHS
    zeta = isolver.BX.dot(zeta)
    zeta = zeta - sum(zeta * w) / sum(w)
    isolver.save_sol[k, :] = zeta
    # zeta = ly.layerpotSD(s=sb, t=so).dot(zeta)
    time.sleep(0.0005)
    isolver.save_ratio[k] = np.sqrt(sum(RHS**2 * w)) / np.sqrt(sum(zeta**2 * w))
  return
def iallsols_one(isolver, w, k, k_alpha):
  # 1st linear sistem
  RHS = isolver.RHS_fcom(isolver.RHS_args)
  RHSr = isolver.RHS_f(isolver, RHS=RHS)
  isolver.Ar = computeAtikh(isolver.A, isolver.alpha)
  zeta = linalg.solve(isolver.Ar, RHSr)
    
  # 2nd linear sistem
  # isolver.Ar = computeAtikh(A, isolver.alpha)
  zeta_p = linalg.solve(isolver.Ar, -zeta)

  res = isolver.A.dot(zeta) - RHS
  res = isolver.BY.dot(res)
  Kdx = isolver.A.dot(zeta_p)
  Kdx = isolver.BY.dot(Kdx)

  # save
  isolver.save_zeta[k, :, k_alpha] = zeta
  sol = isolver.BX.dot(zeta)
  sol = sol - sum(sol * w) / sum(w)
  # zeta = ly.layerpotSD(s=sb, t=so).dot(zeta)
  isolver.save_sol[k, :, k_alpha] = sol
  isolver.save_rhs[k, :, k_alpha] = RHS

  sol_p = isolver.BX.dot(zeta_p)
  sol_p = sol_p - sum(sol_p * w) / sum(w)

  # old
  # disc = sum(res**2 * w)
  # isolver.save_disc[k, k_alpha] = np.sqrt(disc)

  # disc = np.sqrt(disc) - isolver.delta * np.sqrt(sum(zeta**2 * w))
  # disc_p = 2 * sum(res * Kdx * w)
  disc = func_disc_pow1_n(res, sol, w, isolver.delta)
  disc_p = func_disc_p_pow1_n(res, Kdx, sol, sol_p, w, isolver.delta)
  isolver.save_disc[k, k_alpha] = disc

  isolver.save_alpha[k, k_alpha] = isolver.alpha
  isolver.save_disc_p[k, k_alpha] = disc_p
  isolver.save_ratio[k, k_alpha] = np.sqrt(sum(RHS**2 * w)) / np.sqrt(sum(sol**2 * w))
  
  # update alpha
  # isolver.alpha = isolver.alpha - disc / disc_p
  isolver.alpha, isolver.alpha_l, isolver.alpha_r = func_alpha_bis(isolver.alpha, isolver.alpha_l, isolver.alpha_r, disc)
  return
def iallsols_opt(isolver, pointstest, so, it_alpha=2):
  w = so.w

  # ninv = np.empty(len(pointstest.x), float)
  isolver.save_zeta = np.empty((len(pointstest.x), isolver.A.shape[1], it_alpha), float)
  isolver.save_sol = np.empty((len(pointstest.x), so.n, it_alpha), float)
  isolver.save_alpha = np.empty((len(pointstest.x), it_alpha), float)
  isolver.save_disc = np.empty((len(pointstest.x), it_alpha), float)
  isolver.save_disc_p = np.empty((len(pointstest.x), it_alpha), float)
  isolver.save_ratio = np.empty((len(pointstest.x), it_alpha), float)

  # test one point
  isolver.RHS_args['z0'] = pointstest.x[0]
  RHS = isolver.RHS_fcom(isolver.RHS_args)

  isolver.save_rhs = np.empty((len(pointstest.x), len(RHS), it_alpha), float)
  alpha_orig = isolver.alpha
  for k in range(len(pointstest.x)):
    isolver.RHS_args['z0'] = pointstest.x[k]
    isolver.alpha = alpha_orig
    isolver.alpha_l = 1e-16
    isolver.alpha_r = alpha_orig
    for k_alpha in range(it_alpha):
      iallsols_one(isolver, w, k, k_alpha)
      # time.sleep(0.005)
  return
#####################################################
def eigselect(A, m0 = ()):
  As = 0.5 * (A + A.T)
  if max([max(abs(r)) for r in As]) > 1e-10:
    print('Warning: not sym matrix in eigselect')
  (w, v) = linalg.eig(As)
  wind = [(abs(w)[k], k) for k in range(len(w))]
  wind = sorted(wind, key=lambda x: x[0], reverse=True)
  if m0 == ():
    m0 = 15
  if len(w) < m0:
    print('Warning: too much eigs')
  windt = list(map(list, zip(*wind)))
  plt.plot(range(len(w)), np.log(windt[0]), 'b-+')
  x = range(m0)
  y = np.log(windt[0][0:m0])
  linreg = scipy.stats.linregress(x, y)
  plt.plot(x, linreg.intercept + linreg.slope*x, 'r-+')
  plt.show(block=False)
  return(w, v, wind, m0, linreg)

def ieig(w, v, wind, m0, linreg, isolver, pointstest, LL0, so, theta=0):
  windt = list(map(list, zip(*wind)))
  res = ()
  nsolgap = ()
  ninv = np.empty(len(pointstest.x), float)

  weigths = ()
  if len(w) == so.n:
    print('scalar product n')
    weigths = so.w
  elif len(w) == so.n - 1:
    print('scalar product n - 1')
    weigths = np.ones(so.n - 1)
  else:
    print('Some error in length of w')
  for k in range(len(pointstest.x)):
    isolver.RHS_args['z0'] = pointstest.x[k]

    RHS = isolver.RHS_fcom(isolver.RHS_args)
    rhs_coeffs = RHS.T.dot(np.diagflat(weigths).dot(v[:, windt[1][0:m0]]))
    x = range(m0)
    rhs_y = np.array(np.log(abs(rhs_coeffs)**2), float) # needed to avoid Type error in Python (bug?)
    rhs_linreg = scipy.stats.linregress(x, rhs_y)

    time.sleep(0.001)
    if np.exp(rhs_linreg.slope) < np.exp(linreg.slope):
      ninv[k] = 1
    else:
      ninv[k] = 0
    ninv[k] = np.exp(rhs_linreg.slope) - np.exp(linreg.slope)
  return (ninv, res, nsolgap)

####################################
def test_one(isolver, pointstest, so, alpha):
  w = so.w

  ninv = np.empty(len(pointstest.x), float)
  isolver.save_zeta = np.empty((len(pointstest.x), isolver.A.shape[1], len(alpha)), float)
  isolver.save_sol = np.empty((len(pointstest.x), so.n, len(alpha)), float)
  isolver.save_alpha = np.empty((len(pointstest.x), len(alpha)), float)
  isolver.save_disc = np.empty((len(pointstest.x), len(alpha)), float)
  isolver.save_disc_p = np.empty((len(pointstest.x), len(alpha)), float)
  isolver.save_ratio = np.empty((len(pointstest.x), len(alpha)), float)

  # test one point
  isolver.RHS_args['z0'] = pointstest.x[0]
  RHS = isolver.RHS_fcom(isolver.RHS_args)

  isolver.save_rhs = np.empty((len(pointstest.x), len(RHS), len(alpha)), float)
  # alpha = np.array([al_start + al_step * k for k in range(len(alpha))])
  alpha_orig = isolver.alpha
  for k in range(len(pointstest.x)):
    isolver.RHS_args['z0'] = pointstest.x[k]
    isolver.alpha_l = 1e-16
    isolver.alpha_r = alpha_orig

    for k_alpha in range(len(alpha)):
      isolver.alpha = alpha[k_alpha]
      iallsols_one(isolver, w, k, k_alpha)
      time.sleep(0.005)
  return
