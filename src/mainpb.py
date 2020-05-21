import numpy as np
import matplotlib.pyplot as plt
import scipy.linalg as linalg
import numpy.linalg
#-------------------------------
import layerpot as ly
import segment as sg
import shapes as sh
import geometries as gm
import plot
import directproblem as dpb
import inverseproblem as ipb
import mainpb as  m
# import problem as pb
import linfunc as linf
import lintype as lint


import setups
import numpy.random

import debug_globals as dbgg
# tt = time.time()
# tc = time.clock()

regmet = 'tikh'
solver = 'lu'

BY_set = 0

print('Y regularization with method:', regmet, ' and solver:', solver)

class EIT:
  def __init__(self, alpha = 1e-10, delta = 5e-8, theta = 0, noiselevel = 0, noiselevelK = 0, noiselevelU = 0, m0 = 30, m = 15, h = 2.0): # try noiselevel 0.01
    '''
    This constructor initializes these fields
      - meshgrid: only for final intensity visualization
      - alpha: regularization parameter?
      - delta: ?
      - theta: direction of rhs test function (quite arbitrary)
      - noiselevel
      - noiselevelK
      - c: constant computed from conductivity h
      - p = pointset is set of testing point z0 in the grid
    '''
    self.p = ()
    self.alpha = alpha
    self.delta = delta
    self.theta = theta
    self.noiselevel = noiselevel
    self.noiselevelK = noiselevelK
    self.noiselevelU = noiselevelU
    self.m = m
    self.m0 = m0
    self.h = h
    self.c = sg.const_from_conductivity(1.0 * (h + 1) / (h - 1))
  def domain(self, index='one_ellipse', nsb=80, nso=80, nsd=40, rb=5, ro=3, quad='ps', *args, **kwargs):
    '''
    This function add fields 
      'sb', 'so', 'ld'
    to instance using template function with circumferences for 'sb', 'so'
    '''
    self.sb, self.so, self.ld = gm.template_bod(name_sd=index, nsb=nsb, nso=nso, nsd=nsd, rb=rb, ro=ro, quad=quad, *args, **kwargs)
    rounded_range = np.ceil(np.round(np.max(np.abs(self.so.x)), decimals=2))
    print("--- rounded range = %d from %.2f---" %(rounded_range, np.max(np.abs(self.so.x))))
    self.meshgrid_args = (-rounded_range, rounded_range, 80)
    self.meshgrid()
  def meshgrid(self, args=(), args_y=((), (), ())):
    '''
      meshgrid computes x, y, pp(points plot), p(points computation) (pp is p by default)
    '''
    if args == ():
      args = self.meshgrid_args
    else:
      self.meshgrid_args = args
    self.x, self.y, self.testpoints = ipb.meshgrid(args, args_y)
    self.flag_points_inside(self.so)
    self.z = np.zeros(self.testpoints.flag_inside_s.shape)
  def flag_points_inside(self, s):
    for k in range(len(self.testpoints.x)):
      self.testpoints.flag_inside_s[k] = s.contains(self.testpoints.x[k])
  # -------------------------------------------------------------------
  def plot_domain(self):
    self.ld.plot(p=True)
    self.so.plot()
    plt.show(block=False)
  def plot_pre(self, default_value = 0):
    for k in range(len(self.testpoints.x)):
      if self.testpoints.flag_inside_s[k] == 0:
        self.z[k] = default_value
  def plot(self, z=(), pltype='im', pllevels='exp', *args, **kwargs):
    if z == ():
      z = self.z
    self.plot_pre()
    fig = plot.plot(self.x, self.y, z, pltype=pltype, pllevels=pllevels, *args, **kwargs)
    self.plot_domain()
    plt.show(block=False)
  # -------------------------------------------------------------------
  def addnoiseK(self, K, noiselevel):
    noiseamplitude = np.max(np.abs(K)) * noiselevel
    noise = noiseamplitude * numpy.random.normal(0, 1, K.shape)
    K = K + noise
    print('')
    print('  added noise K')
    return K
  # -------------------------------------------------------------------
  def basis_X(self, t='e'):
    self.so.BX, self.so.BXinv = sg.get_Basis(n=self.so.n, w=self.so.w, t=t)
  def basis_Y(self, t='e'):
    self.so.BY, self.so.BYinv = sg.get_Basis(n=self.so.n, w=self.so.w, t=t)
  def basis_XY(self, t='e'):
    self.basis_X(t)
    self.basis_Y(t)
  # =========================================================================================================================================
  # ----------------------------------------------------- linear sampling -------------------------------------------------------------------
  # =========================================================================================================================================
  def solver(self):
    print('')
    print('--------------------------------------------------------')
    print('-------- linear sampling solver method of EIT ----------')
    print(' - alpha (solver) =', self.alpha)
    if BY_set == 0:
      print(' - BY_set: no (solver)')
      self.L0 = ipb.computeL0(self.so, ())
      self.L0B = ()
      self.L = ipb.computeL(self.ld, self.so, (), self.c)
      self.LL0 = ipb.computeLL0(self.ld, self.so, T=self.so.BX, c=self.c, L0=self.L0, L=self.L)
      RHS_fcom = ipb.NtoD_computeRHS
      K = self.LL0
    elif BY_set == 1:
      print(' - BY_set: yes (solver)')
      self.L0 = ()
      self.L0B = ipb.computeL0B(self.so, ())
      self.LB = ipb.computeLB(self.ld, self.so, (), self.c)
      self.LL0B = ipb.computeLL0B(self.ld, self.so, T=self.so.BX, c=self.c, L0B=self.L0B, LB=self.LB)
      RHS_fcom = ipb.NtoD_computeRHSB
      K = self.LL0B
    # self.u, self.w, self.v = linalg.svd(self.K)
    BXr = np.eye(self.so.n)
    BXr = linf.base_mean(BXr, self.so.w)
    L0 = ipb.computeL0(self.so, BXr)
    L = ipb.computeL(self.ld, self.so, BXr, self.c)
    
    RHS_args = {'L0' : self.L0, 'L0B' : self.L0B, 's' : self.so, 'z0' : (), 'theta' : self.theta, 'noiselevel' : self.noiselevel}
    K = self.addnoiseK(K, self.noiselevelK)
    self.K = K
    self.isolver = ipb.solver_init(A=K, alpha=self.alpha, delta=self.delta, reg=1, regmet=regmet, solvertype=solver, RHS_fcom=RHS_fcom, RHS_args=RHS_args, so=self.so, testpoints=self.testpoints, BX=self.so.BX, BY=self.so.BY)
  # =========================================================================================================================================
  #----------------------------------------------------- reciprocity gap --------------------------------------------------------------------
  # =========================================================================================================================================
  def rg_solver(self):
    print('')
    print('--------------------------------------------------------')
    print('------- reciprocal gap solver method of EIT ------------')
    if BY_set == 0:
      print('BY_set: no')
      self.L0 = ipb.computeL0(self.so, ())
      self.L0B = ()
      self.L = ipb.computeL(self.ld, self.so, (), self.c)
      self.LL0 = ipb.computeLL0(self.ld, self.so, T=self.so.BX, c=self.c, L0=self.L0, L=self.L)
      RHS_fcom = ipb.NtoD_computeRHS
      # ----------------------------------------------------------------
      allrhs, allpsi = ipb.computeallpsi(self.ld, self.sb, self.c)
      # dpb.plotdpb(ld, sb.x[0], (-2, 2, 100), psi = allpsi[:,0], t='im')
      R, U, U_nu = ipb.computeR(allpsi, self.ld, self.so, self.sb)
      RHS_fcom = ipb.gap_computeRHS
    elif BY_set == 1:
      print('BY_set: yes')
    # -------------------------------------------------------------------------------------------------------------------
    '''
      Setting:
      - Solver.RHS_args: 'R', 'U', 'U_nu', 'w'
    '''
    U    = self.addnoiseK(U, self.noiselevelU)
    U_nu = self.addnoiseK(U_nu, self.noiselevelU)
    RHS_args = {'L0' : (), 'L0B' : (), 's' : self.so, 'z0' : (), 'theta' : self.theta, 'noiselevel' : self.noiselevel}
    self.isolver = ipb.solver_init(A=R, alpha=self.alpha, delta=self.delta, reg=1, regmet=regmet, solvertype=solver, RHS_fcom=RHS_fcom, RHS_args=RHS_args, so=self.so, testpoints=self.testpoints, BX=(), BY=())
    self.isolver.RHS_args['R'] = R
    self.isolver.RHS_args['U'] = U
    self.isolver.RHS_args['U_nu'] = U_nu
    # self.isolver.w = self.sb.w, fbxbbzzbz
    return
  def solver_update_alpha(self, alpha_new):
    self.isolver.update_alpha_original(alpha_new)
  # =========================================================================================================================================
  # ---------------------------------------------------- heavy methods : morozov iterations -----------------------------------------------
  # =========================================================================================================================================
  def ipb(self):
    ''' it solves iallsols for alpha fixed '''
    ipb.iallsols(self.isolver, self.testpoints, self.so)  
    # ---- save output ----
    self.z = self.isolver.save_ratio[:, 0]
    self.plot_pre()
  def ipb_alpha_fixedarray(self, alpha=()):
    self.isolver.alpha_fixedarray = np.array([10**(-k) for k in range(10, 50, 5)]) if alpha == () else alpha
    self.isolver.alpha_r_orig, self.isolver.alpha_l_orig, self.isolver.alpha_orig = self.isolver.alpha_fixedarray[0], self.isolver.alpha_fixedarray[-1], self.isolver.alpha_fixedarray[0]
    self.isolver.alpha_method = 'fixedarray'
    ipb.iallsols_opt_new(self.isolver, self.testpoints, len(self.isolver.alpha_fixedarray))
    # ---- save output ----
    self.z = self.isolver.save.ratio[-1, :]
    self.plot_pre()
  def ipb_opt(self, it_alpha = 2):
    '''
      it solves iallsols_opt with Morozov principle and save ratio in z
      - fix number of iterations = it_alpha
    '''
    ipb.iallsols_opt_new(self.isolver, self.testpoints, it_alpha)
    # ---- save output ----
    self.z = self.isolver.save.ratio[-1, :]
    self.plot_pre()
  # =========================================================================================================================================
  def visual_new_discrepancy(self, new_delta, showplot=False):
    self.new_delta = new_delta
    res = (self.isolver.save.zeta.dot(self.isolver.A.T)  - self.isolver.save.rhs).dot(self.so.BY.T)
    self.isolver.save.new_disc = ipb.func_disc_value_pow1_n(res, self.isolver.save.sol, self.so.w, self.new_delta)
    # ---- plot ----
    if showplot:
      plt.figure(); plt.title('new delta = %e' %self.new_delta)
      for k, testpoint_inside in enumerate(self.testpoints.flag_inside_s):
        if testpoint_inside:
          plt.plot(self.isolver.save.alpha[:, k], self.isolver.save.new_disc[:, k])
      plt.plot([0, self.isolver.alpha_r_orig], [0, 0], 'k--')
      plt.show(block=False)
  def visual_vanish_new_discrepancy(self):
    idx_sorted = np.argsort(self.isolver.save.alpha, axis=0)
    self.sorted_alpha = np.take_along_axis(self.isolver.save.alpha, idx_sorted, axis=0)
    self.sorted_discr = np.take_along_axis(self.isolver.save.new_disc, idx_sorted, axis=0)
    # ---- check increasing discr ----
    self.plot(np.all((self.sorted_discr[1:, :] - self.sorted_discr[:-1, :]) >= 0, axis=0), vmin=0, vmax=1)
    plt.title('increasing discrepancy'); plt.show(block=False)
    nbins = self.sorted_alpha.shape[1]
    idx_zero = np.argmax(np.concatenate([np.logical_and(self.sorted_discr[1:, :] > 0, self.sorted_discr[:-1, :] < 0), np.ones((1, nbins), dtype=bool)]), axis=0)
    self.idx_zero = idx_zero
    temp_zero_alpha = np.where(np.logical_or((self.sorted_discr[(idx_zero + 1) % self.sorted_alpha.shape[0], np.arange(nbins)] - self.sorted_discr[idx_zero, np.arange(nbins)]) <= 0, idx_zero == self.sorted_alpha.shape[0]), np.nan, self.sorted_discr[(idx_zero + 1) % self.sorted_alpha.shape[0], np.arange(nbins)] - self.sorted_discr[idx_zero, np.arange(nbins)])
    self.zero_alpha = self.sorted_alpha[idx_zero, np.arange(nbins)] + (-self.sorted_discr[idx_zero, np.arange(nbins)]) / temp_zero_alpha * (self.sorted_alpha[(idx_zero + 1) % self.sorted_alpha.shape[0], np.arange(nbins)] - self.sorted_alpha[idx_zero, np.arange(nbins)])
    
  # =========================================================================================================================================
  #----------------------------------------------------- factorization ----------------------------------------------------------------------
  # =========================================================================================================================================
  def fact_solver(self):
    self.solver()
    self.LLdiff = ipb.computeLLdiff(self.ld, self.so, T=np.eye(self.so.n), c=c)
    self.LLBdiff = ipb.computeLLBdiff(self.ld, self.so, T=self.so.BX[:, 1:], c=c)

    RHS_args = {'L0' : self.L0, 'L0B' : self.L0B, 's' : self.so, 'z0' : (), 'theta' : self.theta, 'noiselevel' : self.noiselevel}
    RHS_fcom = ipb.NtoD_computeRHS
    self.isolver = ipb.solver_init(A=self.LLdiff, alpha=self.alpha, delta=self.delta, reg=1, regmet=regmet, solvertype=solver, RHS_fcom=RHS_fcom, RHS_args=RHS_args, so=self.so, testpoints=self.testpoints, BX=self.so.BX, BY=self.so.BY)
    self.LLfact = self.LL0
    # self.LLfact = self.LLdiff
    K = self.LLfact
    K = self.addnoiseK(K, self.noiselevelK)
    self.K = K
    self.LLfact = self.K
  def fact_addnoiseK(self, K):
    K = self.LLfact
    K = self.addnoiseK(K, self.noiselevelK)
    self.K = K
    self.LLfact = self.K
  def fact_ieig(self):
    mselect = self.m0
    if setups.fact_L_trigbasis:
      BT_red = self.so.BT[:, :self.m]
      self.LLfact = BT_red.T.dot(self.LLfact.dot(BT_red))
      mselect = self.m
    w, v, wsorted, m0, linreg = ipb.eigselect(self.LLfact, m0=mselect)
    # ipb.eigplot(wsorted, m0, linreg)
    self.fact_wsorted = wsorted
    self.fact_linreg = linreg
    print('m0000 = ', m0)
    self.z = ipb.ieig(w, v, wsorted, m0, linreg, self.isolver, self.testpoints, self.LL0, self.so, self.theta)
    self.plot_pre()
