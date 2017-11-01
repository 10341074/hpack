from src import *

def mapNtoD(lo, ld, g, c, s0=()):
  no = lo.n
  nd = ld.n
  Kpd = ly.layerpotSD(s=ld)
  Kpo = ly.layerpotSD(s=lo)
  Kpd[np.diag_indices(nd)] = Kpd[np.diag_indices(nd)] + 0.5 * c
  Kpo[np.diag_indices(no)] = Kpo[np.diag_indices(no)] + 0.5
  Kd2o = ly.layerpotSD(s=ld, t=lo)
  Ko2d = ly.layerpotSD(s=lo, t=ld)

  if s0 == ():
    s0 = np.ones(no)
  S = linf.gramschmidt(s0 = s0)
  Kpo = Kpo.dot(S)
  Ko2d = Ko2d.dot(S)
  
  row1 = np.concatenate((Kpo.T, Kd2o.T)).T
  row2 = np.concatenate((Ko2d.T, Kpd.T)).T
  Ks = np.concatenate((lo.SL.T.dot(row1), row2))
  # Ks = np.concatenate(( row1, row2 ))
  Ks1 = Ks[1::, 1::]

  (lu, piv) = linalg.lu_factor(Ks1)
  if g.ndim == 1:
    gs = np.concatenate((lo.SL.T.dot(g), np.zeros(nd)))
    # gs = np.concatenate(( g, np.zeros(nd) ))
    gs1 = gs[1::]
    if verbose:
      print('mapNtoD condition number= ', numpy.linalg.cond(np.array(Ks, float)))
      print('mapNtoD determninant= ', numpy.linalg.det(np.array(Ks, float)))
    phi1 = linalg.lu_solve((lu, piv), gs1)
    if verbose:
      print('residual = ', numpy.linalg.norm(Ks1.dot(phi1) - gs1))
      print('residual2 = ', numpy.linalg.norm(row2[:, 1::].dot(phi1) - gs1[-nd::]))
    phi = np.concatenate(( [0], phi1 ))
  elif g.ndim == 2:
    nt = g.shape[1]
    gs = np.concatenate(( lo.SL.T.dot(g), np.zeros((nd, nt)) ))
    # gs = np.concatenate(( g, np.zeros((nd, nt)) ))
    gs1 = gs[1::]
    # gs2t = gs2.T
    # phi2 = np.empty((nt, no + nd - 1))
    # for k in range(nt):
    #   phi2[k] = linalg.lu_solve((lu, piv), gs2t[k])
    #   time.sleep(0.001)
    # phi2 = phi2.T
    phi1 = linalg.lu_solve((lu, piv), gs1)
    phi = np.concatenate((np.zeros((1, nt)), phi1))  
  else:
    print('Error dimensions for gs1 in mapNtoD')
  return np.concatenate((S.dot(phi[0:no]), phi[no::]))
def computeL(ld, so, T, c):
  if T == ():
    T = so.BX
  print('computing L')
  allpsi = dpb.mapNtoD(so, ld, T, c, so.s0)
  # allpsi = dpb.mapNtoDD0(so, ld, T, c, so.s0)
  # allpsi = dpb.mapNtoDD_correctedinfirst(so, ld, T, c, so.s0)
  # allpsi = dpb.mapNtoDD_left(so, ld, T, c, so.s0)
  
  Lo = ly.layerpotS(s=so)
  Ld = ly.layerpotS(s=ld, t=so)
  L = Lo.dot(allpsi[0:so.n]) + Ld.dot(allpsi[so.n::])
  # means = sum(np.diagflat(so.w).dot(L)) / sum(so.w) # correct? strange sum by rows
  means = np.ones(so.n).dot(np.diagflat(so.w).dot(L)) / sum(so.w)
  L = L - np.array([means for k in range(so.n)])
  return L
def computeLp(ld, so, T, c, p):
  if T == ():
    T = so.BX
  print('computing Lp')
  allpsi = dpb.mapNtoD(so, ld, T, c, so.s0)
  Lo = ly.layerpotS(s=so)
  Ld = ly.layerpotS(s=ld, t=so)
  L = Lo.dot(allpsi[0:so.n]) + Ld.dot(allpsi[so.n::])
  # means = sum(np.diagflat(so.w).dot(L)) / sum(so.w) # correct? strange sum by rows
  means = np.ones(so.n).dot(np.diagflat(so.w).dot(L)) / sum(so.w)
  ############################
  # new
  Lo = ly.layerpotS(s=so, t=p)
  Ld = ly.layerpotS(s=ld, t=p)
  L = Lo.dot(allpsi[0:so.n]) + Ld.dot(allpsi[so.n::])

  L = L - np.array([means for k in range(Lo.shape[0])])
  return L
def computeL0p(so, T, p):
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
def example():
  p = m.EIT()
  p.domain()
  p.meshgrid((-6,6,40))
  p.z = computeLp(p.ld, p.sb, ly.phi_n(z0=7, z=p.sb.x, n=p.so.nx), 1.02, p.p)
  plot.plot(p.x, p.y, p.z)
  p.plot_domain()
  #############################################################################
  # unu1 = computeL0p(p.sb, ly.phi_n(z0=7, z=p.sb.x, n=p.so.nx), p.so)
  unu1 = computeLp(p.ld, p.sb, ly.phi_n(z0=7, z=p.sb.x, n=p.so.nx), 1.02, p.so)
  unu2 = ly.phi(z0=7, z=p.so.x)
  plt.figure()
  print(p.so.w.dot(unu1)); unu1 = unu1 - sum(p.so.w * unu1) / sum(p.so.w)
  print(p.so.w.dot(unu2)); unu2 = unu2 - sum(p.so.w * unu2) / sum(p.so.w)
  plt.plot(unu1,'+-')
  plt.plot(unu2, '+-')
  plt.show(block=False)
  ############################################################################
def v_ex(z):
  return z.real**3 - 3 * z.real * z.imag**2
def v_p_ex(z):
  return (3 * z.real**2 - 3 * z.imag**2) + 1j * ( - 6 * z.real * z.imag)
def check():
  p = m.EIT()
  p.domain()
  p.meshgrid((-6,6,40))
  r1 = np.empty(p.sb.n)
  phi = ly.scalar(v_p_ex(p.so.x), p.so.nx)
  for k, z0 in enumerate(p.sb.x):
    unu2 = ly.phi(z0=z0, z=p.so.x)
    unu2 = unu2 - sum(p.so.w * unu2) / sum(p.so.w)
    Lphi = ipb.computeLL0(p.ld, p.so, phi, 1.02)
    r1[k] = sum(p.so.w * unu2 * Lphi)
    print(z0)
  p.c = 1.02
  print(p.c)
  p.rg_solver()
  #
  allrhs, allpsi = ipb.computeallpsi(p.ld, p.sb, p.c)
  # dpb.plotdpb(ld, sb.x[0], (-2, 2, 100), psi = allpsi[:,0], t='im')
  R, U, U_nu = ipb.computeR(allpsi, p.ld, p.so, p.sb)
  print("diff matrr = ", tls.mat_max(abs(R- p.K)))
  #
  r2 = p.K.dot(phi)
  plt.figure()
  plt.plot(r1,'+-')
  plt.plot(r2,'+-')
  plt.show(block=False)
  
  
if __name__ == "__main__":
  example()
  check()
  inp = input("Press")
