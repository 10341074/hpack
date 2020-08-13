import matplotlib.pyplot as plt
import numpy as np
from numpy.linalg import norm
from numpy.linalg import cond
import scipy.linalg as linalg
import numpy.linalg
import scipy.stats
import numpy.random
import time

import shapes as sh
import layerpot as ly
import segment as sg
import plot
import directproblem as dpb

import linfunc as linf
import lintype as lint
import setups

from compmapproblem import *

import debug_globals as dbgg
from debug import printname

verbose = 0
saving_disc = 1
# -------------------------------------------------------------
def printname_ipb(func):
  return printname(func, "ipb", print_this=verbose)
def printname_ipb_solve_testpoint(func):
  return printname(func, "ipb", color=dbgg.GREEN, print_this=verbose)
# -------------------------------------------------------------
def meshgrid(x1_x2_xn, y1_y2_yn=((), (), ())):
  x, y, pp = plot.meshgrid(x1_x2_xn, y1_y2_yn)
  pp = sg.Pointset(pp)
  return (x, y, pp)
# ================================================================================================
# ---- rhs: reciprocal gap ----
@printname_ipb
def gap_computeRHS(args):
  '''
    reciprocal gap method: rhs
      - U, U_nu: functions u and derivates u_nu for source points on boundary
      - z0: test point inside domain
      - s: segment, usually so (intermediate sd<so<sb) where defined gap functional
      - (theta: phi parameter)
    this function computes with phi
      - F, F_nu: from phi
  '''
  U, U_nu, z0, so, theta = args['U'], args['U_nu'], args['z0'], args['s'], args['theta']
  F = np.cos(theta) * ly.phi_x(z0, so.x) + np.sin(theta) * ly.phi_y(z0, so.x)
  F_nu = np.cos(theta) * ly.phi_x_n(z0, so.x, so.nx) + np.sin(theta) * ly.phi_y_n(z0, so.x, so.nx)
  rhs = U.dot(F_nu) - U_nu.dot(F)
  # -------------------------------------------------------
  # -- adding noise --
  noiselevel     = args['noiselevel']
  noiseamplitude = max(abs(rhs)) * noiselevel
  noise = noiseamplitude * numpy.random.normal(0, 1, rhs.size)
  # ------------------
  rhs = rhs + noise
  # -------------------------------------------------------
  return rhs

@printname_ipb
def gap_computeRHSB(args):
  rhs = gap_computeRHS(args, rhs)
  return s.B[:, 1:].T.dot(np.diagflat(s.w).dot(rhs))
# ================================================================================================
# -------------------------------------------------------------
# ---- Tikhnov method to solve A * x = rhs, with alpha ----
@printname_ipb
def computeRHStikh(A, rhs):
  '''
    rhs_tikh = A.T * rhs
  '''
  return A.T.dot(rhs)
@printname_ipb
def computeAtikh(A, alpha):
  '''
    A_tikh = alpha * (I + A.T * A)
  '''
  return alpha * np.eye(A.shape[1]) + A.T.dot(A)
# ================================================================================================
# ------------------------ initializations --------------------------------------------------------
@printname_ipb
def NtoD_init(LL0, a, reg, regmet, solver):
  if reg or reg == 1:
    if regmet == 'tikh':
      Rreg = computeRtikh(LL0, a)
      _A_b = (Rreg, computeRHSNtoDtikh)
  else:
    _A_b = (LL0, computeRHSNtoD)

  if solver == 's':
    _gap = (lint._solve, _A_b[0], _A_b[1])
  elif solver == 'lu':
    (lu, piv) = linalg.lu_factor(_A_b[0])
    _gap = (lint._lu, (lu, piv), _A_b[1])
  elif solver == 'lstsq':
    _gap = (lint._lstsq, _A_b[0], _A_b[1])
  else:
    print('Error: not valid solver')
  return _gap

@printname_ipb
def solver_init(A, alpha, delta, reg, regmet, solvertype, RHS_fcom, RHS_args, so, testpoints, BX=(), BY=()):
  '''
  Questa funzione inizializza e restituisce Solver object. Richiede:
    - tipo di regolarizzazione: per determinare funzioni che regolarizzano matrice A e rhs
    - BX, BY
    - parametri per calcolare rhs (in futuro): rhs_args, funzione da chiamare per rhs (in futuro)
    - regolarizzazione: alpha
  In futuro:
    - manca solo da calcolare rhs, e regolarizzarlo
  '''
  print(dbgg.GREEN, end='')
  print('')
  print('---------------------------------------------------------------------------')
  print('-------------------------- building Solver --------------------------------')
  print(dbgg.RESET, end='')

  if regmet == 'tikh':
    Areg_func   = computeAtikh
    rhsreg_func = computeRHStikh
  else:
    raise Exception('Error: regularization without specified method')
  
  s = lint.Solver(A=A, Areg_func=Areg_func, RHS_f=rhsreg_func, RHS_fcom=RHS_fcom, RHS_args=RHS_args, BX=BX, BY=BY, alpha=alpha, delta=delta, w=so.w, testpoints=testpoints, solvertype=solvertype)
  return s
# ================================================================================================
# ------------------------ reciprocity gap -------------------------------------------------------
@printname_ipb
def computeU(allpsi, ld, so, sb):
  kerS  = ly.layerpotS(s=ld, t=so)
  kerSD = ly.layerpotSD(s=ld, t=so)

  U_psi    = kerS.dot(allpsi)
  U_psi_nu = kerSD.dot(allpsi)

  U_psi    = U_psi.T.dot(np.diag(so.w))
  U_psi_nu = U_psi_nu.T.dot(np.diag(so.w))

  U    = ly.layerpotS(s=so, t=sb)
  U_nu = ly.layerpotD(s=so, t=sb)
  
  U    = U + U_psi
  U_nu = U_nu + U_psi_nu
  return (U, U_nu)

@printname_ipb
def computeR(allpsi, ld, so, sb, sv=(), testset=1):
  if sv == ():
    sv = sb
  nso, nsb = so.n, sb.n

  U, U_nu = computeU(allpsi, ld, so, sb)
  
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
    print('  R.shape = ', R.shape)
  return (R, U, U_nu)
# ================================================================================================
# ------------------------ giant inverse loop ----------------------------------------------------
@printname_ipb
def computeallpsi(ld, sb, c):
  nsd = ld.n
  nsb = sb.n
  allpsi = np.empty((nsd, nsb), float)
  allrhs = np.empty((nsd, nsb), float)
  A = dpb.directKpc(ld, c)
  for k in range(nsb):
    z0 = sb.x[k]
    allrhs[:, k] = dpb.directrhs(l=ld, z0=z0)
  for k in range(nsb):
    allpsi[:, k] = dpb.directpb(l=ld, c=c, z0=z0, A_f = (A, linalg.solve), rhs=allrhs[:, k])
  return allrhs, allpsi

# ================================================================================================
# ---------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------
@printname_ipb
def computeRHSNtoD(LL0, z0, so, theta=0, RHS=()):
  if RHS != ():
    return RHS
  S = linf.gramschmidtw(so)
  RHS = S[:, 1::].T.dot(ly.phi_theta(z0, so.x, theta))
  return RHS

@printname_ipb
def computeRHSNtoDtikh(LL0, z0, so, theta=0, RHS=()):
  if RHS == ():
    RHS = computeRHSNtoD(LL0, z0, so, theta)
  return LL0.T.dot(RHS)
# ================================================================================================
@printname_ipb
def NtoD_computeneumb(args):
  ''' neumann term in linear sampling rhs '''
  s, z0, theta = args['s'], args['z0'], args['theta']
  xx = ly.phi_xx(z0, s.x)
  xy = ly.phi_xy(z0, s.x)
  yy = ly.phi_yy(z0, s.x)
  d = np.cos(theta) + 1j * np.sin(theta)
  neum = s.nx.real * xx * d.real + s.nx.real * xy * d.imag + \
         s.nx.imag * xy * d.real + s.nx.imag * yy * d.imag
  return neum

# ================================================================================================
# ---- rhs: linear sampling ----
@printname_ipb
def NtoD_computeRHS(args, rhs=()):
  '''
    linear sampling method: rhs = (3 terms) = rhspsi - m (mean) - dirich(to_vanish_norm_deriv)
  '''
  L0, L0B, s, z0, theta = args['L0'], args['L0B'], args['s'], args['z0'], args['theta']
  if rhs != (): # remove check at every call
    return rhs
  ''' term 3 '''
  neumb = NtoD_computeneumb(args)
  dirh = L0.dot(neumb)
  ''' term 1 '''
  d = np.cos(theta) + 1j * np.sin(theta)
  rhs = ly.scalar(d, ly.phi_p(z0, s.x))
  ''' term 2 '''
  m = sum(rhs * s.w) / sum(s.w)
  rhs = rhs - m
  rhs = rhs - dirh
  # -------------------------------------------------------
  # -- adding noise --
  noiselevel     = args['noiselevel']
  noiseamplitude = max(abs(rhs)) * noiselevel
  noise = noiseamplitude * numpy.random.normal(0, 1, rhs.size)
  # ------------------
  rhs = rhs + noise
  # -------------------------------------------------------
  return rhs

@printname_ipb
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
  rhs = rhs - dirh
  # --------------
  return rhsh
# ================================================================================================
def func_disc_value_pow1_n(res, sol, w, delta):
  '''
    The following formula is required to vanish
      | K * x{alpha, delta} - y{delta} | -  delta |x{alpha, delta}|
      the norm is interpeted as an integration
    The parameters are
    - res = residual
    - sol = x
  '''
  ret = np.sqrt((res**2).dot(w)) - delta * np.sqrt((sol**2).dot(w))
  # print("  discrepancy shape =", ret.shape)
  return ret
def func_disc_deriv_pow1_n(res, Kdx, sol, sol_p, w, delta):
  '''
    The following formula (x is shortend for x{delt, alpha}, y for y{delta}) is first variation applied to increment dx
      1/| K*x - y| * 2 * <K*x - y, K*dx)> - delta * 1/|x| * 2 * <x, dx>
  '''  
  ret = 2 * (res * Kdx).dot(w) / np.sqrt((res**2).dot(w)) - delta * 2 * (sol * sol_p).dot(w) / np.sqrt((sol**2).dot(w))
  # print("  discrepancy deriv shape =", ret.shape)
  return ret
# -----------------------------------------------------
def func_alpha_newton(alpha, discrepancy, discrepancy_derived):
  return alpha - discrepancy / discrepancy_derived
def func_alpha_bis(alpha, alpha_l, alpha_r, discrepancy):
  if discrepancy > 0:
    alpha_r = alpha
  else:
    alpha_l = alpha
  alpha = (alpha_l + alpha_r) * 0.5
  return (alpha, alpha_l, alpha_r)
# ================================================================================================
@printname_ipb
def NtoD_computeRHS__tikh__print():
  return
@printname_ipb_solve_testpoint
def iallsols(isolver, pointstest, so):
  '''
    compute all sols, named zeta, of inverse problem: (N_reg) zeta = rhs
  '''
  print(dbgg.CYAN, end='')
  print('  condition num = %.5e' %numpy.linalg.cond(np.array(isolver.A, float)))
  print('  determinant   = %.5e' %numpy.linalg.det(np.array(isolver.A, float)))
  print(dbgg.RESET, end='')


  w = isolver.w # previously was set to sb.w (error?)
  if isolver.w == ():
    w = so.w
  isolver.save_zeta  = np.empty((len(pointstest.x), isolver.A.shape[1]), float)
  isolver.save_sol   = np.empty((len(pointstest.x), isolver.A.shape[1]), float)
  isolver.save_ratio = np.empty((len(pointstest.x),1), float)

  # test one point
  isolver.RHS_args['z0'] = pointstest.x[0]
  RHS = isolver.RHS_fcom(isolver.RHS_args)

  isolver.save_rhs = np.empty((len(pointstest.x), len(RHS)), float)
  for k in range(len(pointstest.x)):
    if pointstest.flag_inside_s[k] == 1:
      isolver.RHS_args['z0'] = pointstest.x[k]
      # --- printing ---
      NtoD_computeRHS__tikh__print()
      # ----------------
      dbgg.print_this = False
      RHS = isolver.RHS_fcom(isolver.RHS_args)
      dbgg.print_this = False
      RHSr = isolver.RHS_f(isolver.A, RHS)
      # ==========================================
      # --- solving ----
      zeta = isolver.solver_f(isolver.solver_a, RHSr)
      # ==========================================
      # --- going back to real-valued function ---
      isolver.save_zeta[k, :] = zeta
      isolver.save_rhs[k] = RHS
      zeta = isolver.BX.dot(zeta)
      zeta = zeta - sum(zeta * w) / sum(w)
      # ---- finished computations ----------------------------------------
      # ---- saving: zeta (real values on so.x), ratio (from zeta, and rhs)
      isolver.save_sol[k, :] = zeta
      isolver.save_ratio[k] = np.sqrt(sum(RHS**2 * w)) / np.sqrt(sum(zeta**2 * w))
  return
@printname_ipb_solve_testpoint
def iallsols_onetestpoint(isolver, w, k_point, k_alpha, rhs):
  '''
    First equation in zeta:
      eq 1.0: (alpha + K^T K) zeta = K^T y
    Second equation in d(zeta)/d(alpha) derives from differentiation of first equation with respect to alpha:
    alpha * zeta + K^T K zeta = K^T y
    zeta + alpha d(zeta)/d(alpha) + K^T K d(zeta)/d(alpha) = 0
      eq 2.0: (alpha + K^T K) d(zeta)/d(alpha) = - zeta
  '''
  # ---- 1st linear sistem -----------------------------------
  rhsregul   = isolver.RHS_f(isolver.A, rhs=rhs)
  isolver.Ar = computeAtikh(isolver.A, isolver.alpha)
  # print("isolver for alpha %e cond %e " %( isolver.alpha, numpy.linalg.cond(np.array(isolver.Ar, float))))
  zeta = linalg.solve(isolver.Ar, rhsregul)

  # ---- 2nd linear sistem -----------------------------------
  zeta_p = linalg.solve(isolver.Ar, -zeta)
  # ----------------------------------------------------------
  # ---- save densities of linear equation -------
  sol = isolver.BX.dot(zeta)
  sol = sol - sum(sol * w) / sum(w)
  if saving_disc:
    isolver.save.zeta[k_alpha, k_point, :] = zeta
    isolver.save.sol [k_alpha, k_point, :] = sol
    isolver.save.rhs [k_alpha, k_point, :] = rhs

  sol_p = isolver.BX.dot(zeta_p)
  sol_p = sol_p - sum(sol_p * w) / sum(w)
  # ----------------------------------------------------------
  res = isolver.BY.dot(isolver.A.dot(zeta) - rhs)
  Kdx = isolver.BY.dot(isolver.A.dot(zeta_p))
  # ---- save ------------------------------------------------
  discrepancy       = func_disc_value_pow1_n(res, sol, w, isolver.delta)
  discrepancy_deriv = func_disc_deriv_pow1_n(res, Kdx, sol, sol_p, w, isolver.delta)
  if saving_disc:
    isolver.save.disc    [k_alpha, k_point] = discrepancy
    isolver.save.alpha   [k_alpha, k_point] = isolver.alpha
    isolver.save.discr_d [k_alpha, k_point] = discrepancy_deriv
    isolver.save.ratio   [k_alpha, k_point] = np.sqrt(sum(rhs**2 * w)) / np.sqrt(sum(sol**2 * w))
    isolver.save.alpha_l          [k_point] = isolver.alpha_l
    isolver.save.alpha_r          [k_point] = isolver.alpha_r
  if k_alpha == (isolver.save.ratio.shape[0] - 1):
    isolver.save.ratio   [k_alpha, k_point] = np.sqrt(sum(rhs**2 * w)) / np.sqrt(sum(sol**2 * w))
  # ----------------------------------------------------------
  # ---- new alpha limits for next iterations -----
  if isolver.alpha_method == 'opt_bisect':
    isolver.alpha, isolver.alpha_l, isolver.alpha_r = func_alpha_bis(isolver.alpha, isolver.alpha_l, isolver.alpha_r, discrepancy)
  elif isolver.alpha_method == 'fixedarray':
    isolver.alpha, isolver.alpha_l, isolver.alpha_r = isolver.alpha_fixedarray[(k_alpha+1)%len(isolver.alpha_fixedarray)], isolver.alpha_fixedarray[-1], isolver.alpha_fixedarray[0]
  else:
    raise Exception('Error: not valid method for alpha')
  return
@printname_ipb_solve_testpoint
def iallsols_opt_new(isolver, pointstest, it_alpha):
  print(dbgg.CYAN, end='')
  print('  matrix Ar with alpha = %e' %isolver.alpha)
  print('  condition num = %.5e' %numpy.linalg.cond(np.array(isolver.Ar, float)))
  print('  determinant   = %.5e' %numpy.linalg.det(np.array(isolver.Ar, float)))
  print(dbgg.RESET, end='')
  isolver.add_iterations_alpha(num_testpoints=len(pointstest.x), it_alpha=it_alpha)
  # -----------------------------------------------------------------------------
  for k_point in range(len(pointstest.x)):
    if pointstest.flag_inside_s[k_point] == 1:
      isolver.RHS_args['z0'] = pointstest.x[k_point]
      isolver.alpha   = isolver.alpha_orig
      isolver.alpha_l = isolver.alpha_l_orig
      isolver.alpha_r = isolver.alpha_r_orig
      rhs = isolver.RHS_fcom(isolver.RHS_args) # to put outside
      if pointstest.flag_inside_s[k_point] == 1:
        for k_alpha in range(it_alpha):
          iallsols_onetestpoint(isolver, isolver.w, k_point, k_alpha, rhs = rhs)
          # time.sleep(0.005)
    if np.mod(k_point, 50) == 0:
      print('  percentage %.2f' %(k_point / len(pointstest.x) * 100))
  return
###########################################################################################################
#### FACTORIZATION
###########################################################################################################
def eigselect(A, m0 = 15):
  # symmetrize matrix
  As = 0.5 * (A + A.T)
  if max([max(abs(r)) for r in As - As.T]) > 1e-10:
    print('WARNING: not sym matrix in eigselect')
  # compute eigenvalues
  (w, v) = linalg.eig(As)
  # sort eigenvalues
  wsortedt = [(abs(w)[k], k) for k in range(len(w))]
  wsortedt = sorted(wsortedt, key=lambda x: x[0], reverse=True)
  if len(w) < m0:
    print('WARNING: not enough eigenvalues')
  # transposition
  wsorted = list(map(list, zip(*wsortedt)))
  # linear regression
  x = np.arange(setups.fact_w_discard, m0 - setups.fact_w_discard)
  yk = np.array(wsorted[0][setups.fact_w_discard:m0 - setups.fact_w_discard])
  if setups.fact_meaned:
    x = x[: int(x.size / 2)]
    # yk = np.concatenate(( [yk[0]], np.sqrt(yk[1:-1:2] * yk[2:-1:2]) ))
    yk = np.concatenate(( [yk[0]], 1.0 / (1.0 / yk[1:-1:2] + 1.0 / yk[2:-1:2]) ))
  y = np.log(yk)
  linreg = scipy.stats.linregress(x, y)
  # plot
  # fig = plt.figure()
  # plt.plot(range(len(w)), np.log(windt[0]), 'b-+')
  # plt.plot(x, linreg.intercept + linreg.slope*x, 'r-+')
  # plt.show(block=False)
  return(w, v, wsorted, m0, linreg)

def eigonly(A, m0 = 15):
  w, v, wsorted, m0, linreg = eigselect(A, m0)
  return wsorted
def eigplot(wsorted, m0, linreg, savefig = 0):
  # plot
  fig = plt.figure()
  plt.plot(range(len(wsorted[0])), np.log(wsorted[0]), 'bp', ms=2)
  x = np.arange(m0)
  plt.plot(x, linreg.intercept + linreg.slope*x, 'r-*', ms=0.2)
  plt.show(block=False)  
  return
def ieig(w, v, wsorted, m0, linreg, isolver, pointstest, LL0, so, theta=0):
  # transposition
  # windt = list(map(list, zip(*wind)))
  # initialization
  chi = np.empty(len(pointstest.x), float)
  # initialize rhs scalar product
  weigths = ()
  if len(w) == so.n:
    print('check: rhs components = n')
    weigths = so.w
  elif len(w) == so.n - 1:
    print('check: rhs componenents = n - 1')
    weigths = np.ones(so.n - 1)
  else:
    print('ERROR: some error in length of w')
  # loop
  for k in range(len(pointstest.x)):
    isolver.RHS_args['z0'] = pointstest.x[k]
    RHS = isolver.RHS_fcom(isolver.RHS_args)
    if setups.fact_L_trigbasis:
      BT_red = so.BT[:, :m0]
      RHS = BT_red.T.dot(RHS) 
    rhs_coeffs = RHS.T.dot(v[:, wsorted[1][0:m0]])
    # rhs_coeffs = RHS.T.dot(np.diagflat(weigths).dot(v[:, wsorted[1][0:m0]]))
    # linear regression
    x = np.arange(setups.fact_w_discard, m0 - setups.fact_w_discard)
    yk = rhs_coeffs[setups.fact_w_discard : m0 - setups.fact_w_discard]
    if setups.fact_meaned:
      x = x[: int(x.size / 2)]
      yk = np.concatenate(( [yk[0]], yk[1:-1:2] + yk[2:-1:2] ))
    rhs_y = np.array(np.log(abs(yk)**2), float) # needed to avoid Type error in Python (bug?)
    rhs_linreg = scipy.stats.linregress(x, rhs_y)

    if rhs_linreg.slope < linreg.slope:
      chi[k] = 1
    else:
      chi[k] = 0
    # print(pointstest.x[k], chi[k])
    # ninv[k] = np.exp(rhs_linreg.slope) - np.exp(linreg.slope)
  return chi
############################################
def testieig(w, v, wsorted, m0, linreg, isolver, pointstest, LL0, so, theta=0):
  # transposition
  # windt = list(map(list, zip(*wind)))
  # initialization
  chi = np.empty(len(pointstest.x), float)
  # initialize rhs scalar product
  weigths = ()
  if len(w) == so.n:
    print('check: rhs components = n')
    weigths = so.w
  elif len(w) == so.n - 1:
    print('check: rhs componenents = n - 1')
    weigths = np.ones(so.n - 1)
  else:
    print('ERROR: some error in length of w')
  # loop
  for k in range(len(pointstest.x)):
    isolver.RHS_args['z0'] = pointstest.x[k]
    RHS = isolver.RHS_fcom(isolver.RHS_args)
    if setups.fact_L_trigbasis:
      BT_red = so.BT[:, :m0]
      RHS = BT_red.T.dot(RHS) 
    rhs_coeffs = RHS.T.dot(v[:, wsorted[1][0:m0]])
    # rhs_coeffs = RHS.T.dot(np.diagflat(weigths).dot(v[:, wsorted[1][0:m0]]))
    # linear regression
    x = np.arange(setups.fact_w_discard, m0 - setups.fact_w_discard)
    yk = rhs_coeffs[setups.fact_w_discard : m0 - setups.fact_w_discard]
    if setups.fact_meaned:
      x = x[: int(x.size / 2)]
      yk = np.concatenate(( [yk[0]], yk[1:-1:2] + yk[2:-1:2] ))
    rhs_y = np.array(np.log(abs(yk)**2), float) # needed to avoid Type error in Python (bug?)
    print(rhs_y)
    rhs_linreg = scipy.stats.linregress(x, rhs_y)

    if rhs_linreg.slope < linreg.slope:
      chi[k] = 1
    else:
      chi[k] = 0
    # ninv[k] = np.exp(rhs_linreg.slope) - np.exp(linreg.slope)
  print(rhs_y)
  return chi, rhs_y, rhs_linreg
# ================================================================================================
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
