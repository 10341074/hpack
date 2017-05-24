from __load__ import *
# from inverseproblem import *
import inverseproblem as ipb
#import __mainpb__ as _m # for test

tt = time.time()
tc = time.clock()
#########

h = 2
nsrc = 20

R = ()
alpha = 1e-10
a = alpha # to be deleted
reg = 1
regmet = 'tikh'
solver = 'lu'

theta =  0

BY_set = 0
#########
c = 1.0 * (h + 1) / (h - 1)
print('c = ', c)
print('theta = ', theta)
if reg or reg == 1:
  print('Y regularization with method:', regmet, ' and solver:', solver)
else:
  print('N regularization with solver ', solver)

def x3domain(nsb=0, nso=0, nsd=0):
  ns = nsrc
  
  rb  = 10
  nsb = rb * ns
  if nsb == 0:
    nsb = 160
  sb = sg.Segment(nsb, f_inargs = (sh.circle, (0, rb)), quad='ps')

  ro  = 3
  nso = ro * ns
  if nso == 0:
    nso = 180
  so = sg.Segment(nso, f_inargs = (sh.circle, (0, ro)), quad='ps')
  # so = sg.Segment(nso, f_inargs = (sh.ellipse, (0, 4, 3)), quad='ps')

  rd  = 1
  if nsd == 0:
    nsd = 100
  # sd = sg.Segment(nsd, f_inargs = (sh.circle, (0, rd)), quad='ps')
  sd = sg.Segment(nsd, f_inargs = (sh.ellipse, (0, 2*rd, rd)), quad='ps')
  # sd = sg.Segment(nsd, Z_args = (sh.dZ, sh.dZp, sh.dZpp, ()), quad='gp', aff=(0, 0.8 + 0.8j))
  # sd = sg.Segment(nsd, Z_args = (sh.kZ, sh.kZp, sh.kZpp, ()), quad='ps', aff=(0, 0.8 + 0.8j))
  sd1 = sg.Segment(50, f_inargs = (sh.ellipse, (0, 2*rd, rd)), quad='ps', aff=(-1-1j, 0.4 + 0.4j))
  sd2 = sg.Segment(50, f_inargs = (sh.ellipse, (0, 2*rd, rd)), quad='ps', aff=(0.5 +1j, 0.4 - 0.4j))
  
  
  bd1 = sg.Boundary([sd1])
  bd2 = sg.Boundary([sd2])
  ld = sg.Layer([bd1, bd2], ly.layerpotSD)
  return (ld, so, sb)

def method_gap():
  ld, so, sb = x3domain()
  nsd, nso, nsb = ld.n, so.n, sb.n

  allpsi = computeallpsi(ld, sb, c)
  # dpb.plotdpb(ld, sb.x[0], (-2, 2, 100), psi = allpsi[:,0], t='im')

  R, U, U_nu = computeR(allpsi, ld, so, sb)
  
  RHS_args = {'R': R , 'U': U, 'U_nu': U_nu, 'z0': (), 'so': so, 'theta': theta}
  RHS_fcom = ipb.gap_computeRHSB
  isolver_gap = ipb.solver_init(R, a, reg, regmet, solver, RHS_fcom=RHS_fcom, RHS_args=RHS_args, BX=so.BX, BY=so.BY)
  x, y, pp = meshgrid((-2, 2 , 40))

  (ninv, sol, rhs, res, nsolgap) = ipb.iallsols(isolver_gap, pp, so=so)
 
  plot.plot(x, y, ninv, 'im')
  ld.plot(p=True)
  plt.show(block=False)

  # R = np.array(R, float)
  # Rreg = np.array(Rreg, float)

  # plt.savefig('fig.pdf')
  # plt.savefig('fig.png')
  # plt.savefig('fig.ps')
  # plt.savefig('fig.eps')

  plt.savefig('fig_ninv.svg')
  return (ninv, res)

def method_NtoD(p=()):
  ld, so, sb = x3domain()
  nsd, nso, nsb = ld.n, so.n, sb.n

  # so.BX = ly.layerpotSDnow(s=sb, t=so)
  # so.BX = np.eye(so.n)
  so.BX = so.B[:, 1:]
  # so.BY = ly.layerpotSD(s=sb, t=so)
  # so.BX = linf.base_mean(so.BX, so.w)
  # so.BX = linf.base_norm(so.BX, so.w)

  # so.BX = so.BX[:, 1:]
  # so.BX = linf.base_mean(so.BX, so.w)
  
  if BY_set == 0:
    print('BY_set: no')
    L0 = ipb.computeL0(so, ())
    L0B = ()
    L = ipb.computeL(ld, so, (), c)
    LL0 = ipb.computeLL0(ld, so, T=so.BX, c=c, L0=L0, L=L)
    RHS_fcom = ipb.NtoD_computeRHS
    M = LL0
  elif BY_set == 1:
    print('BY_set: yes')
    L0 = ()
    L0B = ipb.computeL0B(so, ())
    LB = ipb.computeLB(ld, so, (), c)
    LL0B = ipb.computeLL0B(ld, so, T=so.BX, c=c, L0B=L0B, LB=LB)
    RHS_fcom = ipb.NtoD_computeRHSB
    M = LL0B
  
  BXr = np.eye(so.n)
  BXr = linf.base_mean(BXr, so.w)
  L0 = ipb.computeL0(so, BXr)
  L = ipb.computeL(ld, so, BXr, c)

  RHS_args = {'L0' : L0, 'L0B' : L0B, 's' : so, 'z0' : (), 'theta' : theta} 
  isolver_NtoD = ipb.solver_init(M, a, reg, regmet, solver, RHS_fcom=RHS_fcom, RHS_args=RHS_args, BX=so.BX, BY=so.BY)

  x, y, pp = ipb.meshgrid((-2, 1, 20))
  if p == ():
    p = pp
  else:
    p = sg.Pointset(p)
  ipb.iallsols(isolver_NtoD, p, so)  
  # ipb.iallsols_opt(isolver_NtoD, pp, so, it_alpha=2) 

  # ipb.test_one(isolver_NtoD, pp, so, 10, 1e-16, 1e-16) 
  ninv = isolver_NtoD.save_ratio[:, 0]
  res = ()
  plot.plot(x, y, ninv, 'im')
  ld.plot(p=True)
  so.plot()
  plt.show(block=False)
  plt.savefig('fig_ninvL0.svg')
  return isolver_NtoD
def method_F():
  ld, so, sb = x3domain()
  nsd, nso, nsb = ld.n, so.n, sb.n

  LL0 = ipb.computeLL0(ld, so, T=so.B[:, 1:], c=c)
  # LL0B = ipb.computeLL0B(ld, so, T=so.B[:, 1:], c=c, testset=0)
  L0 = ipb.computeL0(so, ())
  L0B = ipb.computeL0B(so, ())

  LLdiff = ipb.computeLLdiff(ld, so, T=np.eye(so.n), c=c)
  LLBdiff = ipb.computeLLBdiff(ld, so, T=so.B[:, 1:], c=c)

  RHS_args = {'L0' : L0, 'L0B' : L0B, 's' : so, 'z0' : (), 'theta' : theta}
  RHS_fcom = ipb.NtoD_computeRHS
  isolver_NtoD = ipb.solver_init(LLdiff, a, reg, regmet, solver, RHS_fcom=RHS_fcom, RHS_args=RHS_args, BX=so.BX, BY=so.BY)
  
  x, y, pp = ipb.meshgrid((-2, 2, 40))
  w, v, wind, m0, linreg = ipb. eigselect(LLdiff, m0=40)
  
  (ninv, res, nsolgap) = ipb.ieig(w, v, wind, m0, linreg, isolver_NtoD, pp, LL0, so, theta)
 
  plot.plot(x, y, ninv, 'im')
  ld.plot(p=True)
  so.plot()
  plt.show(block=False)
  # plt.savefig('fig_ninv.svg')
  return (ninv, res)
#############
ld, so, sb = x3domain()
nsd, nso, nsb = ld.n, so.n, sb.n
if __name__ == "__main__":
  isolver = method_F()
thetav = np.pi * np.array([0], float)
for theta in np.array(thetav):
  # ninv = method_NtoD()
  pass
tt = time.time() - tt
tc = time.clock() - tc

print('time wall-clock = ', tt)
print('time clock = ', tc)

end = input('Press enter')
############
def test_NtoD():
  g = ly.phi_n(-8, so.x, so.n)
  #g = np.ones(so.n)
  psi0 = dpb.mapNtoD0(so, g)
  bo = sg.Boundary([so])
  lo = sg.Layer([bo], dns=psi0)
  
  dpb.plotdpb(lo, (), (-8, 8, 100), t='im')

  psi = dpb.mapNtoD(so, ld, g, c)
  lo.dns = psi[0:lo.n]
  ld.dns = psi[lo.n:]

  dpb.plotdpb(lo, (), (-8, 8, 100), t='im', l2=ld)
  return

def resid(sol, rhs, A):
  r = np.empty(len(sol), float)
  for k in range(len(sol)):
    r[k] = linalg.norm(A.dot(sol[k]) - rhs[k])
  return r

def f(sb, so, ld):
  if BY_set == 0:
    print('BY_set: no')
    L0 = ipb.computeL0(so, so.B[:, 1:])
    L0B = ()
    L = ipb.computeL(ld, so, so.B[:, 1:], c)
    LL0 = ipb.computeLL0(ld, so, T=so.BX, c=c, testset=0)
    RHS_fcom = ipb.NtoD_computeRHS
    M = LL0
  elif BY_set == 1:
    print('BY_set: yes')
    L0 = ()
    L0B = ipb.computeL0B(so, ())
    LB = ipb.computeLB(ld, so, (), c)
    LL0B = ipb.computeLL0B(ld, so, T=so.BX, c=c, testset=0)
    RHS_fcom = ipb.NtoD_computeRHSB
    M = LL0B

  RHS_args = {'L0' : L0, 'L0B' : L0B, 's' : so, 'z0' : (), 'theta' : theta} 
  isolver_NtoD = ipb.solver_init(M, a, reg, regmet, solver, RHS_fcom=RHS_fcom, RHS_args=RHS_args)

  x, y, pp = ipb.meshgrid((-2, 2, 40))
  (ninv, sol, rhs, res, nsolgap) = ipb.iallsols(isolver_NtoD, pp, sb, so)
 
  plot.plot(x, y, ninv, 'cf')
  ld.plot(p=True)
  so.plot()
  plt.show(block=False)
  plt.savefig('fig_ninv.svg')
  return (sol, rhs)
  
def method_test():
  ld, so, sb = x3domain()
  nsd, nso, nsb = ld.n, so.n, sb.n

  # so.BX = ly.layerpotSDnow(s=sb, t=so)
  # so.BX = np.eye(so.n)
  # so.BX = so.B[:, 1:]
  # so.BY = ly.layerpotSD(s=sb, t=so)
  # so.BX = linf.base_mean(so.BX, so.w)
  # so.BX = linf.base_norm(so.BX, so.w)

  so.BX = so.BX[:, 1:]
  so.BX = linf.base_mean(so.BX, so.w)
  
  if BY_set == 0:
    print('BY_set: no')
    L0 = ipb.computeL0(so, ())
    L0B = ()
    L = ipb.computeL(ld, so, (), c)
    LL0 = ipb.computeLL0(ld, so, T=so.BX, c=c, L0=L0, L=L)
    RHS_fcom = ipb.NtoD_computeRHS
    M = LL0
  elif BY_set == 1:
    print('BY_set: yes')
    L0 = ()
    L0B = ipb.computeL0B(so, ())
    LB = ipb.computeLB(ld, so, (), c)
    LL0B = ipb.computeLL0B(ld, so, T=so.BX, c=c, L0B=L0B, LB=LB)
    RHS_fcom = ipb.NtoD_computeRHSB
    M = LL0B
  
  BXr = np.eye(so.n)
  BXr = linf.base_mean(BXr, so.w)
  L0 = ipb.computeL0(so, BXr)

  RHS_args = {'L0' : L0, 'L0B' : L0B, 's' : so, 'z0' : (), 'theta' : theta} 
  isolver_NtoD = ipb.solver_init(M, a, reg, regmet, solver, RHS_fcom=RHS_fcom, RHS_args=RHS_args, BX=so.BX, BY=so.BY)

  x, y, pp = ipb.meshgrid((-2, 1, 10))
  # (ninv, sol, rhs, res, nsolgap) = ipb.iallsols(isolver_NtoD, pp, sb=sb, so=so)
  
  # ipb.iallsols_opt(isolver_NtoD, pp, so, it_alpha=3) 
  ipb.test_one(isolver_NtoD, pp, so, 20, 1e-16, 5) 
  ninv = isolver_NtoD.save_ratio[:, 1]
  res = ()
  plot.plot(x, y, ninv, 'im')
  ld.plot(p=True)
  so.plot()
  plt.show(block=False)
  # return (ninv, res)
  return isolver_NtoD

class EIT:
  def __init__(self, a = a, theta = theta):
    self.meshgrid_args = (-2, 1, 20)
    self.p = ()
    self.alpha = alpha
    self.theta = theta
    return
  def domain(self, nsb=0, nso=0, nsd=0):
    self.ld, self.so, self.sb = x3domain(nsb, nso, nsd)
  def meshgrid(self, args=(), args_y=((), (), ())):
    if args == ():
      args = self.meshgrid_args
    else:
      self.meshgrid_args = args
    self.x, self.y, self.pp = ipb.meshgrid(args, args_y)
    if self.p == ():
      self.p = self.pp
  def plot(self, z=(), t='im'):
    if z == ():
      z = self.z
    plot.plot(self.x, self.y, z, t)
    self.ld.plot(p=True)
    self.so.plot()
    plt.show(block=False)
  def plot_domain(self):
    self.ld.plot(p=True)
    self.so.plot()
    plt.show(block=False)
  def solver(self):
    self.alpha = alpha
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
    
    RHS_args = {'L0' : self.L0, 'L0B' : self.L0B, 's' : self.so, 'z0' : (), 'theta' : self.theta} 
    self.isolver = ipb.solver_init(self.K, self.alpha, reg, regmet, solver, RHS_fcom=RHS_fcom, RHS_args=RHS_args, BX=self.so.BX, BY=self.so.BY)
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
  def alpha_fixed_ratio(self, k=0):
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
def alpha_fixed_ratio(k=0, s=()):
  zeta = s.isolver.save_zeta[:, :, k]
  zeta = zeta.T
  w = s.so.w
  zeta = s.isolver.BX.dot(zeta)
  # zeta = zeta - sum(zeta * w) / sum(w)
  zeta = zeta - np.array([sum(np.diagflat(w).dot(zeta)) for r in range(zeta.shape[0])])
  # zeta = ly.layerpotSD(s=sb, t=so).dot(zeta)
  RHS = s.isolver.save_rhs[:,:,0].T
  # self.z  = np.sqrt(sum(RHS**2 * w)) / np.sqrt(sum(zeta**2 * w))
  s.z  = np.sqrt(sum(np.diagflat(w).dot(RHS**2))) / np.sqrt(sum(np.diagflat(w).dot(zeta**2)))
  
