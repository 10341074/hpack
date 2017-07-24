import matplotlib.pyplot as plt
from src import *
tt = time.time()
tc = time.clock()
#########

h = 2
nsrc = 20

reg = 1
regmet = 'tikh'
solver = 'lu'

theta =  0

BY_set = 0

pause = 0
#########
c = 1.0 * (h + 1) / (h - 1)
print('c = ', c)
print('theta = ', theta)
if reg or reg == 1:
  print('Y regularization with method:', regmet, ' and solver:', solver)
else:
  print('N regularization with solver ', solver)


tt = time.time() - tt
tc = time.clock() - tc

print('time wall-clock = ', tt)
print('time clock = ', tc)

if __name__ == "__main__":
  end = input('Press enter')

class Problem:
  def __init__(self, alpha = 1e-10, delta = 1e-6, theta = 0, noiselevel = 0, m0 = 40): # try noiselevel 0.01
    self.meshgrid_args = (-2, 2, 20)
    self.p = ()
    self.alpha = alpha
    self.delta = delta
    self.theta = theta
    self.noiselevel = noiselevel
    self.m0 = m0
    return
  def domain(self, index='one_ellipse', nsb=80, nso=80, nsd=40):
    self.sb, self.so, self.ld = gm.example(index=index, nsb=nsb, nso=nso, nsd=nsd)
  def meshgrid(self, args=(), args_y=((), (), ())):
    # meshgrid computes x, y, pp(points plot), p(points computation) (pp is p by default)
    if args == ():
      args = self.meshgrid_args
    else:
      self.meshgrid_args = args
    self.x, self.y, self.pp = ipb.meshgrid(args, args_y)
    self.importp()
  def importp(self, s=()):
    self.p = self.pp
    if s == ():
      s = self.so
    for k in range(len(self.p.x)):
      self.p.flag_inside_s[k] = s.contains(self.p.x[k])
  def plot_pre(self, defualt_value = 0):
    for k in range(len(self.p.x)):
      if self.p.flag_inside_s[k] == 0:
        self.z[k] = default_value
  def plot(self, z=(), t='im'):
    if z == ():
      z = self.z
    self.plot_pre()
    fig = plot.plot(self.x, self.y, z, t)
    self.plot_domain()
    plt.show(block=False)
  def plot_domain(self):
    self.ld.plot(p=True)
    self.so.plot()
    plt.show(block=False)
  def solver(self):
    # self.alpha = alpha
    print('alpha for solver', self.alpha)
    if BY_set == 0:
      print('BY_set: no')
      self.L0 = ipb.computeL0(self.so, ())
      self.L0B = ()
      self.L = ipb.computeL(self.ld, self.so, (), c)
      self.LL0 = ipb.computeLL0(self.ld, self.so, T=self.so.BX, c=c, L0=self.L0, L=self.L)
      RHS_fcom = ipb.NtoD_computeRHS
      self.K = self.LL0
    elif BY_set == 1:
      print('BY_set: yes')
      self.L0 = ()
      self.L0B = ipb.computeL0B(self.so, ())
      self.LB = ipb.computeLB(self.ld, self.so, (), c)
      self.LL0B = ipb.computeLL0B(self.ld, self.so, T=self.so.BX, c=c, L0B=self.L0B, LB=self.LB)
      RHS_fcom = ipb.NtoD_computeRHSB
      self.K = self.LL0B
    self.u, self.w, self.v = linalg.svd(self.K)
    
    BXr = np.eye(self.so.n)
    BXr = linf.base_mean(BXr, self.so.w)
    L0 = ipb.computeL0(self.so, BXr)
    L = ipb.computeL(self.ld, self.so, BXr, c)
    
    RHS_args = {'L0' : self.L0, 'L0B' : self.L0B, 's' : self.so, 'z0' : (), 'theta' : self.theta, 'noiselevel' : self.noiselevel} 
    self.isolver = ipb.solver_init(self.K, self.alpha, self.delta, reg, regmet, solver, RHS_fcom=RHS_fcom, RHS_args=RHS_args, BX=self.so.BX, BY=self.so.BY)
    # self.isolver.alpha = self.alpha
    # self.isolver.alpha_orig = self.isolver.alpha
    # self.isolver.alpha_l_orig = 1e-16
    # self.isolver.alpha_l = self.isolver.alpha_l_orig
  def ipb(self):
    ipb.iallsols(self.isolver, self.p, self.so)  
    # ipb.iallsols_opt(isolver_NtoD, pp, so, it_alpha=2)
    alpha = np.array([1e-16 + 1e-16 * k for k in range(10)])    
    # ipb.test_one(isolver_NtoD, pp, so, 10, 1e-16, 1e-16) 
    self.z = self.isolver.save_ratio[:, 0]
  def test_alpha(self, alpha=()):
    if alpha == ():
      self.alpha = np.array([1e-16 + 1e-16 * k for k in range(10)])
    else:
      self.alpha = alpha
    ipb.test_one(self.isolver, self.p, so, self.alpha) 
    self.z = self.isolver.save_ratio[:, 0]
  def ipb_opt(self, it_alpha = 2):
    print('alpha_orig', self.isolver.alpha)
    ipb.iallsols_opt(self.isolver, self.p, self.so, it_alpha) 
    # ipb.iallsols(self.isolver, self.p, self.so)  
    # alpha = np.array([1e-16 + 1e-16 * k for k in range(10)])    
    # ipb.test_one(isolver_NtoD, pp, so, 10, 1e-16, 1e-16) 
    self.z = self.isolver.save_ratio[:, 0]
  def ipb_opt_append(self, it_alpha = 2):
    print('alpha_orig', self.isolver.alpha)
    ipb.iallsols_opt_append(self.isolver, self.p, self.so, it_alpha) 
    self.z = self.isolver.save_ratio[:, 0]
  def alpha_fixed_ratio(self, k=0):
    # compute ratio for iteration on alpha fixed equal to k
    zeta = self.isolver.save_zeta[:, :, k]
    zeta = zeta.T
    w = self.so.w
    zeta = self.isolver.BX.dot(zeta)
    # zeta = zeta - sum(zeta * w) / sum(w)
    zeta = zeta - np.array([sum(np.diagflat(w).dot(zeta)) for r in range(zeta.shape[0])])
    # zeta = ly.layerpotSD(s=sb, t=so).dot(zeta)
    RHS = self.isolver.save_rhs[:,:,0].T
    # self.z  = np.sqrt(sum(RHS**2 * w)) / np.sqrt(sum(zeta**2 * w))
    self.z  = np.sqrt(sum(np.diagflat(w).dot(RHS**2))) / np.sqrt(sum(np.diagflat(w).dot(zeta**2)))
  def basis_X(self, t='e'):
    self.so.BX, self.so.BXinv = sg.get_Basis(n=self.so.n, w=self.so.w, t=t)
  def basis_Y(self, t='e'):
    self.so.BY, self.so.BYinv = sg.get_Basis(n=self.so.n, w=self.so.w, t=t)
  def basis_XY(self, t='e'):
    self.basis_X(t)
    self.basis_Y(t)
  def fact_solver(self):
    self.solver()
    self.LLdiff = ipb.computeLLdiff(self.ld, self.so, T=np.eye(self.so.n), c=c)
    self.LLBdiff = ipb.computeLLBdiff(self.ld, self.so, T=self.so.B[:, 1:], c=c)
    RHS_args = {'L0' : self.L0, 'L0B' : self.L0B, 's' : self.so, 'z0' : (), 'theta' : self.theta, 'noiselevel' : self.noiselevel}
    RHS_fcom = ipb.NtoD_computeRHS
    self.isolver = ipb.solver_init(self.LLdiff, self.alpha, self.delta, reg, regmet, solver, RHS_fcom=RHS_fcom, RHS_args=RHS_args, BX=self.so.BX, BY=self.so.BY)
  def fact_ieig(self):
    w, v, wind, m0, linreg = ipb.eigselect(self.LLdiff, m0=self.m0)  
    self.z, self.res, nsolgap = ipb.ieig(w, v, wind, m0, linreg, self.isolver, self.pp, self.LL0, self.so, self.theta)

# def alpha_fixed_ratio(k=0, s=()):
#   zeta = s.isolver.save_zeta[:, :, k]
#   zeta = zeta.T
#   w = s.so.w
#   zeta = s.isolver.BX.dot(zeta)
#   # zeta = zeta - sum(zeta * w) / sum(w)
#   zeta = zeta - np.array([sum(np.diagflat(w).dot(zeta)) for r in range(zeta.shape[0])])
#   # zeta = ly.layerpotSD(s=sb, t=so).dot(zeta)
#   RHS = s.isolver.save_rhs[:,:,0].T
#   # self.z  = np.sqrt(sum(RHS**2 * w)) / np.sqrt(sum(zeta**2 * w))
#   s.z  = np.sqrt(sum(np.diagflat(w).dot(RHS**2))) / np.sqrt(sum(np.diagflat(w).dot(zeta**2)))
  
