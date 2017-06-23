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
p = m.EIT()
p.domain()
ro = 3


p.so = sg.Segment(300, f_inargs = (sh.circle, (0, ro)), quad='p')
allpsi0 = dpb.mapNtoD0(p.so, p.so.BX[0], p.so.s0)
a = allpsi0
err =[]
for n in range(100, 301, 4):
  p.so = sg.Segment(n, f_inargs = (sh.circle, (0, ro)), quad='p')
  b = dpb.mapNtoD0(p.so, p.so.BX[0], p.so.s0)
  err.append([b[0] - a[0], b[int(n/4)] - a[50], b[int(n/2)] - a[100], b[int(n*3/4)] - a[150]])
p.so = sg.Segment(200, f_inargs = (sh.circle, (0, ro)), quad='p')
b = dpb.mapNtoD0(p.so, p.so.BX[0], p.so.s0)
# err.append([b[0] - a[0], b[50] - a[50], b[100] - a[100], b[int(n*3/4)] - a[150]])
