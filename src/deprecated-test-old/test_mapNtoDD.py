from src import *

c = 3
def v_ex(z):
  return z.real**3 - 3 * z.real * z.imag**2
def v_p_ex(z):
  return (3 * z.real**2 - 3 * z.imag**2) + 1j * ( - 6 * z.real * z.imag)

def computeL(ld, so, a):
  # allpsi = dpb.mapNtoD(so, ld, T, c, so.s0)
  # allpsi = dpb.mapNtoDD0(so, ld, T, c, so.s0)
  # allpsi = dpb.mapNtoDD_correctedinfirst(so, ld, T, c, so.s0)
  allpsi = a  
  Lo = ly.layerpotS(s=so)
  Ld = ly.layerpotS(s=ld, t=so)
  L = Lo.dot(allpsi[0:so.n]) + Ld.dot(allpsi[so.n::])
  # means = sum(np.diagflat(so.w).dot(L)) / sum(so.w) # correct? strange sum by rows
  means = np.ones(so.n).dot(np.diagflat(so.w).dot(L)) / sum(so.w)
  L = L - np.array([means for k in range(so.n)])
  return L

p = m.EIT()
p.domain('one_kite', nso = 100, nsd = 80)
p.solver()
a1 = dpb.mapNtoD(p.so, p.ld, ly.scalar(v_p_ex(p.so.x), p.so.nx), c, p.so.s0)
d1 = computeL(p.ld, p.so, a1)
a2 = dpb.mapNtoDD_correctedinfirst(p.so, p.ld, ly.scalar(v_p_ex(p.so.x), p.so.nx), c, p.so.s0)
d2 = computeL(p.ld, p.so, a2)
a3 = dpb.mapNtoDD_left(p.so, p.ld, ly.scalar(v_p_ex(p.so.x), p.so.nx), c, p.so.s0)
d3 = computeL(p.ld, p.so, a3)
plt.plot(d1,'+-')
plt.plot(d2, '+-')
plt.plot(d3, '+-')
plt.plot(v_ex(p.so.x), '+-')
plt.show(block=False)
aa1 = dpb.mapNtoD(p.so, p.ld, np.eye(p.so.n), c, p.so.s0)
dd1 = computeL(p.ld, p.so, aa1)
aa2 = dpb.mapNtoDD_correctedinfirst(p.so, p.ld, np.eye(p.so.n), c, p.so.s0)
dd2 = computeL(p.ld, p.so, aa2)
aa3 = dpb.mapNtoDD_left(p.so, p.ld, np.eye(p.so.n), c, p.so.s0)
dd3 = computeL(p.ld, p.so, aa3)



############################################################################
h = 2
c = 1.0 * (h + 1) / (h - 1)
lo, ld = p.so, p.ld
no = lo.n
nd = ld.n
Kpd = ly.layerpotSD(s=ld)
Kpo = ly.layerpotSD(s=lo)
Kpd[np.diag_indices(nd)] = Kpd[np.diag_indices(nd)] + 0.5 * c
Kpo[np.diag_indices(no)] = Kpo[np.diag_indices(no)] + 0.5
Kd2o = ly.layerpotSD(s=ld, t=lo)
Ko2d = ly.layerpotSD(s=lo, t=ld)

row1 = np.concatenate((Kpo.T, Kd2o.T)).T
row2 = np.concatenate((Ko2d.T, Kpd.T)).T
Ks = np.concatenate((row1, row2))
