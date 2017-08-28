# coding: utf-8
n = 100
from __load__ import *
s = sg.Segment(n, f_inargs = (sh.ellipse, (0, 3, 2)), quad='ps')
Kpa = ly.layerpotSD(s=s)
Kpa[np.diag_indices(n)] = Kpa[np.diag_indices(n)] + 0.5
nu, q = linf.eigmaxpowerw(ly.layerpotSD(s=s), s=s)

g = Kpa.dot(np.ones(n))
# check image of Kpa
sum(g * s.w)

# p = linalg.solve(Kpa, g)
# linalg.norm(Kpa.dot(p) - g) # not big error but warning
# (lu, piv) = linalg.lu_factor(Kpa)
# p = linalg.lu_solve((lu, piv), g)
# linalg.norm(Kpa.dot(p) - g) # not big error but warning

S = linf.gramschmidt(s0=q)
# Kpa.dot(S)[:, 0] = 0 # even if orthogonality is euclidean
Kps = S.T.dot(Kpa.dot(S))
Kps_ = Kps[1::, 1::]

gs = S.T.dot(g)
gs_ = gs[1::]

phi_ = linalg.solve(Kps_, gs_)
#(lu, piv) = linalg.lu_factor(Kps2)
# phi2 = linalg.lu_solve((lu, piv), gs2)
phi = np.concatenate(([0], phi_))
# # phi = linalg.solve(Kps, gs) # check error
# phi = (phi, s_j)
phi = S.dot(phi)

u1 = phi
linalg.norm(Kpa.dot(phi) - g) # ok


# 2
S = linf.gramschmidtw(s=s, s0=q)
# Kpa.dot(S)[:, 0] = 0 # even if orthogonality is euclidean
Kps = S.T.dot(np.diagflat(s.w).dot( Kpa.dot(S) ))
Kps_ = Kps[1::, 1::]

gs = S.T.dot(np.diagflat(s.w).dot(g))
gs_ = gs[1::]

phi_ = linalg.solve(Kps_, gs_)
#(lu, piv) = linalg.lu_factor(Kps2)
# phi2 = linalg.lu_solve((lu, piv), gs2)
phi = np.concatenate(([0], phi_))
# # phi = linalg.solve(Kps, gs) # check error
phi = S.dot(phi)

u2 = phi
linalg.norm(Kpa.dot(phi) - g) # ok

# 3
S = linf.gramschmidt(s0=q)
# Kpa.dot(S)[:, 0] = 0 # even if orthogonality is euclidean
Kps = Kpa.dot(S)
Kps_ = Kps[1::, 1::]

gs = g
gs_ = gs[1::]

phi_ = linalg.solve(Kps_, gs_)
#(lu, piv) = linalg.lu_factor(Kps2)
# phi2 = linalg.lu_solve((lu, piv), gs2)
phi = np.concatenate(([0], phi_))
# # phi = linalg.solve(Kps, gs) # check error
# phi = (phi, s_j)
phi = S.dot(phi)


u3 = phi
linalg.norm(Kpa.dot(u3) - g) # ok

