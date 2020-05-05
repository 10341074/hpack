# coding: utf-8
from __load__ import *
x, y, pp = ipb.meshgrid((-4, 4, 40))
n = 40
s = sg.Segment(n, f_inargs = (sh.ellipse, (0, 3, 2)), quad='ps')
g = ly.phi_x(1, s.x)
K = ly.layerpotD(s=s)
K[np.diag_indices(s.n)] = K[np.diag_indices(s.n)]-0.5
linalg.solve(K, -g)
psi = linalg.solve(K, -g)
Kpp = ly.layerpotD(s=s, t = pp)
vv = Kpp.dot(psi)
vv
vv = Kpp.dot(psi)+ly.phi_x(1,pp.x)
plot.plot(x,y,vv)
plot.plot(x,y,vv,'cf')
