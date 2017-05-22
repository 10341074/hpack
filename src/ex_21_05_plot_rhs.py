from __load__ import *
import mainpb as m

def computeL0(so, T, p):
  if T == ():
    T = so.BX
  print('computing L0')
  allpsi0 = dpb.mapNtoD0(so, T, so.s0) 
  Lo = ly.layerpotS(s=so)
  L0 = Lo.dot(allpsi0)
  means = np.ones(so.n).dot(np.diagflat(so.w).dot(L0)) / sum(so.w)
  # new
  Lo = ly.layerpotS(s=so, t=p)
  L0 = Lo.dot(allpsi0)
  # new
  L0 = L0 - np.array([means for k in range(L0.shape[0])])
  return L0

def computeL(ld, so, T, c, p):
  if T == ():
    T = so.BX
  print('computing L')
  allpsi = dpb.mapNtoD(so, ld, T, c, so.s0)
  Lo = ly.layerpotS(s=so)
  Ld = ly.layerpotS(s=ld, t=so)
  L = Lo.dot(allpsi[0:so.n]) + Ld.dot(allpsi[so.n::])
  # means = sum(np.diagflat(so.w).dot(L)) / sum(so.w) # correct? strange sum by rows
  means = np.ones(so.n).dot(np.diagflat(so.w).dot(L)) / sum(so.w)
  # new
  Lo = ly.layerpotS(s=so, t=p)
  Ld = ly.layerpotS(s=ld, t=p)
  L = Lo.dot(allpsi[0:so.n]) + Ld.dot(allpsi[so.n::])

  L = L - np.array([means for k in range(L0.shape[0])])
  return L
def NtoD_computeRHS(args, p, rhs=()):
  L0, L0B, s, z0, theta = args['L0'], args['L0B'], args['s'], args['z0'], args['theta']
  if rhs != (): # remove check at every call
    return rhs
  a = np.cos(theta) + 1j * np.sin(theta)
  # neumb = NtoD_computeneumb(args)
  # print('mean', sum(neumb*s.w))
  # need decomposition BX of neumb !!!!!!!!!!!!!!!!!!!
  # dirh = L0.dot(neumb)
  dirh = 0
  rhs = ly.scalar(a, ly.phi_p(z0, s.x))
  m = sum(rhs * s.w) / sum(s.w)
  # new
  rhs = ly.scalar(a, ly.phi_p(z0, p.x))
  rhs = rhs - m
  return rhs - dirh

p = m.EIT()
p.domain()
p.meshgrid((-2, 2, 40), (-1, 1, 20))
p.solver()
q = [1.]
p.p = sg.Pointset(q)
# p.plot_domain()
p.ipb()

L0 = computeL0(p.so, p.so.BX, p.pp)
z = L0.dot(p.isolver.save_sol[0])
p.plot(z)


RHS_args = {'L0' : L0, 'L0B' : (), 's' : p.so, 'z0' : q, 'theta' : p.theta}
z = NtoD_computeRHS(RHS_args, p.pp, rhs=())
p.plot(z)
