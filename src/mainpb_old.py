# def method_gap():
#   ld, so, sb = x3domain()
#   nsd, nso, nsb = ld.n, so.n, sb.n

#   allpsi = computeallpsi(ld, sb, c)
#   # dpb.plotdpb(ld, sb.x[0], (-2, 2, 100), psi = allpsi[:,0], t='im')

#   R, U, U_nu = computeR(allpsi, ld, so, sb)
  
#   RHS_args = {'R': R , 'U': U, 'U_nu': U_nu, 'z0': (), 'so': so, 'theta': theta}
#   RHS_fcom = ipb.gap_computeRHSB
#   isolver_gap = ipb.solver_init(R, a, delta, reg, regmet, solver, RHS_fcom=RHS_fcom, RHS_args=RHS_args, BX=so.BX, BY=so.BY)
#   x, y, pp = meshgrid((-2, 2 , 40))

#   (ninv, sol, rhs, res, nsolgap) = ipb.iallsols(isolver_gap, pp, so=so)
 
#   plot.plot(x, y, ninv, 'im')
#   ld.plot(p=True)
#   plt.show(block=False)

#   # R = np.array(R, float)
#   # Rreg = np.array(Rreg, float)

#   # plt.savefig('fig.pdf')
#   # plt.savefig('fig.png')
#   # plt.savefig('fig.ps')
#   # plt.savefig('fig.eps')

#   plt.savefig('fig_ninv.svg')
#   return (ninv, res)

# def method_NtoD(p=()):
#   ld, so, sb = x3domain()
#   nsd, nso, nsb = ld.n, so.n, sb.n

#   # so.BX = ly.layerpotSDnow(s=sb, t=so)
#   # so.BX = np.eye(so.n)
#   so.BX = so.B[:, 1:]
#   # so.BY = ly.layerpotSD(s=sb, t=so)
#   # so.BX = linf.base_mean(so.BX, so.w)
#   # so.BX = linf.base_norm(so.BX, so.w)

#   # so.BX = so.BX[:, 1:]
#   # so.BX = linf.base_mean(so.BX, so.w)
  
#   if BY_set == 0:
#     print('BY_set: no')
#     L0 = ipb.computeL0(so, ())
#     L0B = ()
#     L = ipb.computeL(ld, so, (), c)
#     LL0 = ipb.computeLL0(ld, so, T=so.BX, c=c, L0=L0, L=L)
#     RHS_fcom = ipb.NtoD_computeRHS
#     M = LL0
#   elif BY_set == 1:
#     print('BY_set: yes')
#     L0 = ()
#     L0B = ipb.computeL0B(so, ())
#     LB = ipb.computeLB(ld, so, (), c)
#     LL0B = ipb.computeLL0B(ld, so, T=so.BX, c=c, L0B=L0B, LB=LB)
#     RHS_fcom = ipb.NtoD_computeRHSB
#     M = LL0B
  
#   BXr = np.eye(so.n)
#   BXr = linf.base_mean(BXr, so.w)
#   L0 = ipb.computeL0(so, BXr)
#   L = ipb.computeL(ld, so, BXr, c)

#   RHS_args = {'L0' : L0, 'L0B' : L0B, 's' : so, 'z0' : (), 'theta' : theta} 
#   isolver_NtoD = ipb.solver_init(M, a, delta, reg, regmet, solver, RHS_fcom=RHS_fcom, RHS_args=RHS_args, BX=so.BX, BY=so.BY)

#   x, y, pp = ipb.meshgrid((-2, 1, 20))
#   if p == ():
#     p = pp
#   else:
#     p = sg.Pointset(p)
#   ipb.iallsols(isolver_NtoD, p, so)  
#   # ipb.iallsols_opt(isolver_NtoD, pp, so, it_alpha=2) 

#   # ipb.test_one(isolver_NtoD, pp, so, 10, 1e-16, 1e-16) 
#   ninv = isolver_NtoD.save_ratio[:, 0]
#   res = ()
#   plot.plot(x, y, ninv, 'im')
#   ld.plot(p=True)
#   so.plot()
#   plt.show(block=False)
#   plt.savefig('fig_ninvL0.svg')
#   return isolver_NtoD
# def method_F():
#   ld, so, sb = x3domain()
#   nsd, nso, nsb = ld.n, so.n, sb.n

#   LL0 = ipb.computeLL0(ld, so, T=so.B[:, 1:], c=c)
#   # LL0B = ipb.computeLL0B(ld, so, T=so.B[:, 1:], c=c, testset=0)
#   L0 = ipb.computeL0(so, ())
#   L0B = ipb.computeL0B(so, ())

#   LLdiff = ipb.computeLLdiff(ld, so, T=np.eye(so.n), c=c)
#   LLBdiff = ipb.computeLLBdiff(ld, so, T=so.B[:, 1:], c=c)

#   RHS_args = {'L0' : L0, 'L0B' : L0B, 's' : so, 'z0' : (), 'theta' : theta}
#   RHS_fcom = ipb.NtoD_computeRHS
#   isolver_NtoD = ipb.solver_init(LLdiff, a, delta, reg, regmet, solver, RHS_fcom=RHS_fcom, RHS_args=RHS_args, BX=so.BX, BY=so.BY)
  
#   x, y, pp = ipb.meshgrid((-2, 2, 40))

#   w, v, wind, m0, linreg = ipb.eigselect(LLdiff, m0=40)
  
#   (ninv, res, nsolgap) = ipb.ieig(w, v, wind, m0, linreg, isolver_NtoD, pp, LL0, so, theta)
 
#   plot.plot(x, y, ninv, 'im')
#   ld.plot(p=True)
#   so.plot()
#   plt.show(block=False)
#   # plt.savefig('fig_ninv.svg')
#   return (ninv, res)
#############
# ld, so, sb = x3domain()
# nsd, nso, nsb = ld.n, so.n, sb.n
# if __name__ == "__main__":
#   isolver = method_F()
# thetav = np.pi * np.array([0], float)
# for theta in np.array(thetav):
#   # ninv = method_NtoD()
#   pass
# tt = time.time() - tt
# tc = time.clock() - tc

# print('time wall-clock = ', tt)
# print('time clock = ', tc)


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
  isolver_NtoD = ipb.solver_init(M, a, delta, reg, regmet, solver, RHS_fcom=RHS_fcom, RHS_args=RHS_args)

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
  isolver_NtoD = ipb.solver_init(M, a, delta, reg, regmet, solver, RHS_fcom=RHS_fcom, RHS_args=RHS_args, BX=so.BX, BY=so.BY)

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
