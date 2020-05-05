# coding: utf-8
from __load__ import *
ipb.rhs
s = sg.Segment(n, f_inargs = (sh.ellipse, (0, 3, 2)), quad='ps')
n=40
s = sg.Segment(n, f_inargs = (sh.ellipse, (0, 3, 2)), quad='ps')
help(ipb.rhs)
ipb.rhs(so = s, z0 = 1, theta=0)
reload(ipb)
ipb.rhs(1, s, 0)
ipb.rhs(1, s, np.pi/4)
r1 = ipb.rhs(1, s, np.pi/4)
sum(r1)
sum(r1*s.w)
reload(ipb)
r1 = ipb.rhs(1, s, np.pi/4)
sum(r1*s.w)
(np.cos(np.pi/4) + np.sin(np.4)) * -1
(np.cos(np.pi/4) + np.sin(np.pi/4)) * -1
help(ly)
help(ly.phi_x)
r = ly.phi_x(1, s.x)
sum(r * s.w)
help(ly)
s.nx
help(ly)
s.nx.rela
s.nx.real
reload(ipb)
r1, g = ipb.rhs(1, s, np.pi/4)
g
sum(g * s.w)
dpb.mapNtoD0(-g)
help(dpb.mapNtoD0)
dpb.mapNtoD0(l=s, g=g, s0=linf.eigmaxpowerw(s)[1] )
dpb.mapNtoD0(l=s, g=g, s0=linf.eigmaxpowerw(ly.layerpotSD(s=s), s)[1] )
psi = dpb.mapNtoD0(l=s, g=g, s0=linf.eigmaxpowerw(ly.layerpotSD(s=s), s)[1] )
psi = dpb.mapNtoD0(l=s, g=-g, s0=linf.eigmaxpowerw(ly.layerpotSD(s=s), s)[1] )
