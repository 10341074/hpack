from src import *
#######################################################################
# def v_ex(z):
#   return z.real**3 - 3 * z.real * z.imag**2
# def v_p_ex(z):
#   return (3 * z.real**2 - 3 * z.imag**2) + 1j * ( - 6 * z.real * z.imag)
#######################################################################
def v_ex(z):
  return (z.real - 2)**2 - z.imag**2
def v_p_ex(z):
  return 2 * (z.real - 2) - 2j * z.imag
######################################################################
def computeL0(so, T, p):
  if T == ():
    T = so.BX
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
  allpsi = dpb.mapNtoD(so, ld, T, c, so.s0)
  Lo = ly.layerpotS(s=so)
  Ld = ly.layerpotS(s=ld, t=so)
  L = Lo.dot(allpsi[0:so.n]) + Ld.dot(allpsi[so.n::])
  means = np.ones(so.n).dot(np.diagflat(so.w).dot(L)) / sum(so.w)
  # new
  Lo = ly.layerpotS(s=so, t=p)
  Ld = ly.layerpotS(s=ld, t=p)
  L = Lo.dot(allpsi[0:so.n]) + Ld.dot(allpsi[so.n::])
  # new
  L = L - np.array([means for k in range(L.shape[0])])
  return L
#######################################################################
def NtoD_computeneumb(args):
  s, z0, theta = args['s'], args['z0'], args['theta']
  xx = ly.phi_xx(z0, s.x)
  xy = ly.phi_xy(z0, s.x)
  yy = ly.phi_yy(z0, s.x)
  # hess = [[xx, xy] , [xy, yy]]
  d = np.cos(theta) + 1j * np.sin(theta)
  neum = s.nx.real * xx * d.real + s.nx.real * xy * d.imag + \
         s.nx.imag * xy * d.real + s.nx.imag * yy * d.imag
  return neum

def NtoD_computeRHS(args, p, rhs=()):
  L0, L0B, s, z0, theta = args['L0'], args['L0B'], args['s'], args['z0'], args['theta']
  if rhs != (): # remove check at every call
    return rhs
  d = np.cos(theta) + 1j * np.sin(theta)
  neumb = NtoD_computeneumb(args)
  ext_dirh = L0.dot(neumb)
  rhs = ly.scalar(d, ly.phi_p(z0, s.x))
  m = sum(rhs * s.w) / sum(s.w)
  ext_rhs = ly.scalar(d, ly.phi_p(z0, p.x))
  ext_rhs = ext_rhs - m
  ext_rhs = ext_rhs - ext_dirh
  # noiselevel = args['noiselevel']
  # noisemodul = max(abs(rhs)) * noiselevel
  # noise = noisemodul * numpy.random.normal(0, 1, rhs.size)
  # rhs = rhs + noise
  return ext_rhs
##############################################################################
def plot(p, *args):
  p.plot_pre()
  figf = plt.contourf(p.x, p.y, p.z.reshape(p.y.size, p.x.size), *args)
  fig = plt.contour(p.x, p.y, p.z.reshape(p.y.size, p.x.size), *args, colors='k', linewidths=0.3)
  for c in fig.collections:
    c.set_linestyle('solid')
  p.so.plot(ms = 0.8, lw=2)
  plt.axis('square')
  plt.show()
  return

if __name__ == "__main__":
  c = 2
  p = m.EIT()
  p.domain(nsb=150, nso=200, nsd=150)
  p.meshgrid((-3, 3, 200))
  p.solver()
  
  p.z = computeL0(p.so, ly.scalar(v_p_ex(p.so.x), p.so.nx), p.p)
  plot(p, 50)
  p.z = computeL(p.ld, p.so, ly.scalar(v_p_ex(p.so.x), p.so.nx), c, p.p)
  # p.z = v_ex(p.pp.x)
  plot(p, 50)

  z0 = 1j + 1
  L0 = computeL0(p.so, np.eye(p.so.n), p.p)
  RHS_args = {'L0' : L0, 'L0B' : (), 's' : p.so, 'z0' : z0, 'theta' : p.theta}
  p.z = NtoD_computeRHS(RHS_args, p.p)

  
  V = np.linspace(1.0/50, 2, 100)
  V = np.concatenate((list(reversed(-V)), V))
  
  plot(p, V)
  
  L0 = computeL(p.ld, p.so, np.eye(p.so.n), c, p.p)
  RHS_args = {'L0' : L0, 'L0B' : (), 's' : p.so, 'z0' : z0, 'theta' : p.theta}
  p.z = NtoD_computeRHS(RHS_args, p.p)
  plot(p, V)
