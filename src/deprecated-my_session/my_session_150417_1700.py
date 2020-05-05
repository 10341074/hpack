# coding: utf-8
from __load__ import *
n = 100
s = sg.Segment(n, f_inargs = (sh.ellipse, (0, 3, 2)), quad='ps')
# s = sg.Segment(n, f_inargs = (sh.circle, (0, 2)), quad='ps')
E = np.diagflat(np.ones(n))
# L0 = ipb.computeallpsiL0(s, E)
# w = linalg.eigvals(L0)
# plt.plot(w.real, w.imag, 'b*')
# plt.show(block=True)

Kpa = ly.layerpotSD(s=s)
Kpa[np.diag_indices(l.n)] = Kpa[np.diag_indices(l.n)] + 0.5

# nu, psi = linf.eigmaxpowerw(Kpa, s, nit=200)

Sb = linf.gramschmidtBoundary(so = s, s0 = psi)
L0b = ipb.computeallpsiL0(s, Sb[:, 1:])
w = linalg.eigvals(L0b)
plt.plot(w.real, w.imag, '*b')
plt.show(block=True)


l = s
g = Sb[:,0]
Kp = ly.layerpotSD(s=l)

if verbose:
  print('mapNtoD condition number= ', numpy.linalg.cond(np.array(Kpa, float)))
  print('mapNtoD determninant= ', numpy.linalg.det(np.array(Kpa, float)))
  n = len(Kpa)
  psi = np.ones(n)
  S = gramschmidtBoundary(l, psi)
  Kps = S.T.dot(Kpa.dot(S))
  Kps2 = Kps[1::, 1::]

gs = S.T.dot(g)
gs2 = gs[1::]

# phi2 = linalg.solve(Kps2, gs2)
(lu, piv) = linalg.lu_factor(Kps2)
if gs2.ndim == 1:
  phi2 = linalg.lu_solve((lu, piv), gs2)
  phi = np.concatenate(([0], phi2 ))
  if verbose:
    print('residual = ', numpy.linalg.norm(Kps2.dot(phi2) - gs2))
elif gs2.ndim == 2:
  gs2t = gs2.T
  phi2 = np.empty((len(gs2t), n - 1 ))
  for k in range(len(gs2t)):
    phi2[k] = linalg.lu_solve((lu, piv), gs2t[k])
    time.sleep(0.001)
    phi2 = phi2.T
    phi = np.concatenate((np.zeros((1, len(gs2t))), phi2))
else:
  print('Error dimensions for gs2 in mapNtoD0')
  # phi2 = scipy.sparse.linalg.cg(Kps2, gs2)[0]

# phi = linalg.solve(Kps, gs) # check error
return S.dot(phi)
