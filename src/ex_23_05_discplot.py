from __load__ import *
import mainpb as m
import dill

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
p.domain('two_ellipse')
p.meshgrid((-3, 3, 40))
# p.so.BX = linf.base_mean(p.so.BX, p.so.w)
p.solver()
q = [1.]
p.p = sg.Pointset(q)
p.p = p.pp
# # p.plot_domain()
# p.ipb()
alpha = np.array([1e-16 + 1e-16 * k for k in range(10)])
alpha = np.concatenate((alpha, [1e-15*10**k for k in range(8)]))

#p.test_alpha(alpha=alpha)



# L0 = computeL0(p.so, p.so.BX, p.pp)
# z = L0.dot(p.isolver.save_sol[0])
# p.plot(z)

#
p.solver()
it_alpha = 8
p.ipb_opt(it_alpha=it_alpha)
p.alpha_fixed_ratio(0)
p.plot()
#plt.savefig('prova_0.svg')
p.alpha_fixed_ratio(-1)
p.plot()
#plt.savefig('prova_end.svg')

p.plot()
# dill.dump_session('dill.pkl')
# dill.load_session('dill.pkl')

end = input('Press enter')

# # do this first:
# ! sudo apt-get install cpulimit

# from os import getpid
# from resource import setrlimit, RLIMIT_RSS, RLIM_INFINITY, getrusage, RUSAGE_SELF

# # limit CPU: use only 1% of 1 CPU
# pid = getpid()
# ! cpulimit -b -p $pid -c 1 -l 1

# # limit memory
# setrlimit(RLIMIT_RSS, (50*(1024**2),RLIM_INFINITY))
# print getrusage(RUSAGE_SELF).ru_maxrss

# # this will take a lot of time
# % time sum(xrange(10**8))

# # this will fall over with MemoryError
# sum(range(50*100*(1024**2)))
