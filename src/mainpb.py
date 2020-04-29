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
  def __init__(self, alpha = 1e-10, delta = 1e-6, theta = 0, noiselevel = 0, noiselevelK = 0, m0 = 30, m = 15, h = 2.0): # try noiselevel 0.01
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
    self.m = m
    self.m0 = m0
    self.c = 1.0 * (h + 1) / (h - 1)
  def domain(self, index='one_ellipse', nsb=80, nso=80, nsd=40, *args, **kwargs):
    '''
    This function add fields 
      'sb', 'so', 'ld'
    to instance using template function with circumferences for 'sb', 'so'
    '''
    self.sb, self.so, self.ld = gm.bod_template_c(name_sd=index, nsb=nsb, nso=nso, nsd=nsd, *args, **kwargs)
    rounded_range = np.ceil(np.max(np.abs(self.so.x)))
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
    self.x, self.y, self.pp = ipb.meshgrid(args, args_y)
    self.flag_points_inside(self.so)
  def flag_points_inside(self, s):
    self.p = self.pp
    for k in range(len(self.p.x)):
      self.p.flag_inside_s[k] = s.contains(self.p.x[k])
  # -------------------------------------------------------------------
  def plot_pre(self, default_value = 0):
    for k in range(len(self.p.x)):
      if self.p.flag_inside_s[k] == 0:
        self.z[k] = default_value
  def plot_domain(self):
    self.ld.plot(p=True)
    self.so.plot()
    plt.show(block=False)
  def plot(self, z=(), t='im'):
    if z == ():
      z = self.z
    self.plot_pre()
    fig = plot.plot(self.x, self.y, z, t)
    self.plot_domain()
    plt.show(block=False)
  # -------------------------------------------------------------------
  def addnoiseK(self):
    noisemodul = max([max(abs(r)) for r in self.K]) * self.noiselevelK
    noise = noisemodul * numpy.random.normal(0, 1, self.K.shape)
    self.K = self.K + noise
    print('Added noise K')
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
    self.addnoiseK()
    self.isolver = ipb.solver_init(A=K, alpha=self.alpha, delta=self.delta, reg=1, regmet=regmet, solvertype=solvertype, RHS_fcom=RHS_fcom, RHS_args=RHS_args, so=self.so, testpoints=self.pp, BX=self.so.BX, BY=self.so.BY)
  ###############################################################
  def alpha_fixed_ratio(self, k=0):
    # change output: compute ratio for iteration (on alpha) fixed equal to k
    zeta = self.isolver.save_zeta[:, :, k]
    zeta = zeta.T
    w = self.so.w
    zeta = self.isolver.BX.dot(zeta)
    zeta = zeta - np.array([sum(np.diagflat(w).dot(zeta)) for r in range(zeta.shape[0])])
    RHS = self.isolver.save_rhs[:,:,0].T
    # self.z  = np.sqrt(sum(RHS**2 * w)) / np.sqrt(sum(zeta**2 * w))
    # self.z  = np.sqrt(sum(np.diagflat(w).dot(RHS**2))) / np.sqrt(sum(np.diagflat(w).dot(zeta**2)))
    for k in range(self.p.x.size):
      if self.p.flag_inside_s[k] == 1:
        self.z[k] = np.sqrt(sum(w * RHS[:, k]**2)) / np.sqrt(sum(w * zeta[:, k]**2))
      else:
        self.z[k] = 0
  def basis_X(self, t='e'):
    self.so.BX, self.so.BXinv = sg.get_Basis(n=self.so.n, w=self.so.w, t=t)
  def basis_Y(self, t='e'):
    self.so.BY, self.so.BYinv = sg.get_Basis(n=self.so.n, w=self.so.w, t=t)
  def basis_XY(self, t='e'):
    self.basis_X(t)
    self.basis_Y(t)
  ###############################################################
  # def fact_solver(self):
  #   self.solver()    
  #   self.LLdiff = ipb.computeLLdiff(self.ld, self.so, T=np.eye(self.so.n), c=c)
  #   self.LLBdiff = ipb.computeLLBdiff(self.ld, self.so, T=self.so.B[:, 1:], c=c)
  #   RHS_args = {'L0' : self.L0, 'L0B' : self.L0B, 's' : self.so, 'z0' : (), 'theta' : self.theta, 'noiselevel' : self.noiselevel}
  #   RHS_fcom = ipb.NtoD_computeRHS
  #   self.isolver = ipb.solver_init(self.LLdiff, self.alpha, self.delta, reg, regmet, solver, RHS_fcom=RHS_fcom, RHS_args=RHS_args, BX=self.so.BX, BY=self.so.BY)
  ###############################################################
  def fact_solver(self):
    self.solver()
    self.LLdiff = ipb.computeLLdiff(self.ld, self.so, T=np.eye(self.so.n), c=c)
    self.LLBdiff = ipb.computeLLBdiff(self.ld, self.so, T=self.so.BX[:, 1:], c=c)

    RHS_args = {'L0' : self.L0, 'L0B' : self.L0B, 's' : self.so, 'z0' : (), 'theta' : self.theta, 'noiselevel' : self.noiselevel}
    RHS_fcom = ipb.NtoD_computeRHS
    self.isolver = ipb.solver_init(A=self.LLdiff, alpha=self.alpha, delta=self.delta, reg=1, regmet=regmet, solvertype=solver, RHS_fcom=RHS_fcom, RHS_args=RHS_args, so=self.so, testpoints=self.pp, BX=self.so.BX, BY=self.so.BY)
    self.LLfact = self.LL0
    # self.LLfact = self.LLdiff
    self.K = self.LLfact
    self.addnoiseK()
    self.LLfact = self.K
  def fact_addnoiseK(self):
    self.K = self.LLfact
    self.addnoiseK()
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
    self.z = ipb.ieig(w, v, wsorted, m0, linreg, self.isolver, self.pp, self.LL0, self.so, self.theta)
    self.plot_pre()
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
    RHS_args = {'L0' : (), 'L0B' : (), 's' : self.so, 'z0' : (), 'theta' : self.theta, 'noiselevel' : self.noiselevel}
    self.isolver = ipb.solver_init(A=R, alpha=self.alpha, delta=self.delta, reg=1, regmet=regmet, solvertype=solver, RHS_fcom=RHS_fcom, RHS_args=RHS_args, so=self.so, testpoints=self.pp, BX=(), BY=())
    self.isolver.RHS_args['R'] = R
    self.isolver.RHS_args['U'] = U
    self.isolver.RHS_args['U_nu'] = U_nu
    self.isolver.w = self.sb.w
    return
  # =========================================================================================================================================
  # ---------------------------------------------------- heavy methods : morozov iterations -----------------------------------------------
  # =========================================================================================================================================
  def ipb(self):
    ''' it solves iallsols for alpha fixed '''
    ipb.iallsols(self.isolver, self.p, self.so)  
    # ---- save output ----
    self.z = self.isolver.save_ratio[:, 0]
    self.plot_pre()
  def test_alpha(self, alpha=()):
    if alpha == ():
      self.alpha = np.array([1e-16 + 1e-16 * k for k in range(10)])
    else:
      self.alpha = alpha
    ipb.test_one(self.isolver, self.p, so, self.alpha)
    # ---- save output ----
    self.z = self.isolver.save_ratio[:, 0]
    self.plot_pre()
  def ipb_opt(self, it_alpha = 2):
    '''
      it solves iallsols_opt with Morozov principle and save ratio in z
      - fix number of iterations = it_alpha
    '''
    ipb.iallsols_opt_new(self.isolver, self.p, it_alpha)
    # ---- save output ----
    self.z = self.isolver.save.ratio[0, :]
    self.plot_pre()
  def ipb_opt_append(self, it_alpha = 2):
    '''
      - it appends number of iterations = it_alpha to 'ipb_opt'
    '''
    ipb.iallsols_opt_append(self.isolver, self.p, self.so, it_alpha)
    # ---- save output ----
    self.z = self.isolver.save_ratio[:, 0]
    self.plot_pre()
  # =========================================================================================================================================
  def new_discrepancy(self, new_delta):
    self.isolver.save_new_disc = np.empty(self.isolver.save_disc.shape)
    self.new_delta = new_delta
    for i in range(self.isolver.save_disc.shape[0]):
      for j in range(self.isolver.save_disc.shape[1]):
        self.isolver.save_new_disc[i, j] = ipb.func_disc_value_pow1_n(ipb.computeAtikh(self.isolver.A, self.isolver.save_alpha[i, j]).dot(self.isolver.save_zeta[i, :, j]), self.isolver.save_sol[i, :, j], self.so.w, self.new_delta)
