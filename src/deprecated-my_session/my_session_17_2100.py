# coding: utf-8
from __mainpb__ import *
so
so.plot()
plt.show()
n = 40
s = sg.Segment(n, f_inargs = (sh.ellipse, (0, 3, 2)), quad='ps')
reloa(sg)
reload(sg)
s = sg.Segment(n, f_inargs = (sh.ellipse, (0, 3, 2)), quad='ps')
s.K
help(linf)
nu, s0 = linf.eigmaxpowerw()
K = ly.layerpotSD(s=s)
nu, s0 = linf.eigmaxpowerw(K, s)
nu, s00 = linf.eigmaxpower(K)
s0
s00
s0 / s00
help(linf)
reload(sg)
s = sg.Segment(n, f_inargs = (sh.ellipse, (0, 3, 2)), quad='ps')
reload(sg)
s = sg.Segment(n, f_inargs = (sh.ellipse, (0, 3, 2)), quad='ps')
s.S
s.S[:,0]
s.Sw[:,0]
s.Sw[:,0] / s.S[:, 0]
reload(ipb)
L0 = ipb.computeallpsiL0(s, s.Sw[:, 1::])
L0.shape
A = L0.dot(s.Sw[:, 1::].T.dot(np.diagflat(s.w)))
A
L0 = ipb.computeallpsiL0(s, np.diagflat(np.ones(n)))
n
L0
issym(L0)
sum(L0.T)
sum(np.diagflat(so.w).dot(L0).T)
sum(np.diagflat(so.w).dot(L0))
sum(np.diagflat(s.w).dot(L0).T)
L0.shape
reload(ipb)
L0, m = ipb.computeallpsiL0(s, np.diagflat(np.ones(n)))
L0.shape
reload(ipb)
L0, m = ipb.computeallpsiL0(s, np.diagflat(np.ones(n)))
m
reload(ipb)
L0, m = ipb.computeallpsiL0(s, np.diagflat(np.ones(n)))
sum(np.diagflat(s.w).dot(L0).T)
E = np.diagflat(np.ones(n))
E
allpsi0 = dpb.mapNtoD0(s, E, s.s0)
Lo = ly.layerpotS(s=s)
Lo = Lo.dot(allpsi0)
sum(Lo[:,0]*s.w)
Lo
sum(Lo[:,1]*s.w)
Lo[:, 0]
reload(dpb)
L0 = ipb.computeallpsiL0(s, np.diagflat(np.ones(n)))
reload(dpb)
L0 = ipb.computeallpsiL0(s, np.diagflat(np.ones(n)))
reload(dpb)
E
s.n
M = 1. / np.sum(s.w)
M
M = 1. / np.sum(s.w) * np.array([s.w for k in range(n)])
M
printmat(M)
printmat(M,[1,2],[1,2])
printmat(M,[0::],[1,2])
printmat(M,range(n),[1,2])
printmat(M,range(n),[1,2, 3, 4])
printmat(M,range(n),[1,2, 3])
P = E - M
P.dot(E[:,0])
sum(P.dot(E[:,0]) * s.w)
Lo = ly.layerpotS(s=s)
reload(ipb)
L, m = ipb.computeallpsiL0(s, E)
L
Lo
Lo = ly.layerpotS(s=s)
Lo
allpsi0
dpb.mapNtoD0(s,E[:,0], s.s0)
E[:,0]
ly.layerpoSD(s=s).dot(np.ones(n))
ly.layerpotSD(s=s).dot(np.ones(n))
sum(ly.layerpotSD(s=s).dot(np.ones(n)) * s.w)
sum((ly.layerpotSD(s=s) + 0.5*np.diagflat(np.ones(n))).dot(np.ones(n)) * s.w)
(ly.layerpotSD(s=s) + 0.5*np.diagflat(np.ones(n))).dot(np.ones(n))
q = (ly.layerpotSD(s=s) + 0.5*np.diagflat(np.ones(n))).dot(np.ones(n))
q
sol1 = dpb.mapNtoD0(s,q, s.s0)
alls = dpb.mapNtoD0(s, E, s.s0)
alls.dot(q)
sol2 = alls.dot(q)
sol2
sol2 - sum(sol2*s.w)/sum(s.w)
sol2
sol2 -1
sol1 - sum(sol1*s.w)/sum(s.w)
sol2 - sum(sol2*s.w)/sum(s.w)
ss1=sol1 - sum(sol1*s.w)/sum(s.w) 
ss2 = sol2 - sum(sol2*s.w)/sum(s.w)
ss1 -ss2
alls = dpb.mapNtoD0(s, E, s.s0)
Lo = ly.layerpotS(s=s)
L = Lo.dot(alls)
L.shape
means = np.ones((1,n)).dot(np.diagflat(s.w).dot(L))
means
sum(np.diagflat(s.w).dot(L).T)
A=np.array([[1,2],[0,0]])
sum(A)
A
means
sum(np.diagflat(s.w).dot(L))
reloa(ipb)
reload(ipb)
L = ipb.computeallpsiL0(s,E)
reload(ipb)
L = ipb.computeallpsiL0(s,E)
L
issym(L)
L
L = ipb.computeallpsiL0(s,E)[0]
issym(L)
L
issym(L)
L = ipb.computeallpsiL0(s,E)[0]
L.shape
L,m = ipb.computeallpsiL0(s,E)
L.shape
reload(ipb)
L,m = ipb.computeallpsiL0(s,E)
m
reload(ipb)
L,m = ipb.computeallpsiL0(s,E)
m
L.shape
issym(L)
issym(np.diagflat(s.w).dot(L))
L.shape
sol1 = dpb.mapNtoD0(s,q, s.s0)
sol1
sol1 - sum(sol1*s.w)/sum(s.w)
L.dot(q)
sum(np.diagflat(s.w).dot(L)[:,1])
sum(np.diagflat(s.w).dot(L)[:,2])
sum(np.diagflat(s.w).dot(L)[:,3])
q
sum(np.diagflat(s.w).dot(q))
ly.layerS(s=s).dot(dpb.mapNtoD0(s,q, s.s0))
ly.layerpotS(s=s).dot(dpb.mapNtoD0(s,q, s.s0))
i=ly.layerpotS(s=s).dot(dpb.mapNtoD0(s,q, s.s0))
i - sum(i*s.w)/sum(s.w)
L.dot(q)
L
a = dpb.mapNtoD0(s, E[:,0],s.s0)
a
q
K
K = ly.layerpotSD(s=s) + 0.5*np.eye(n)
np.eye(n)
s0
K.dot(s0)
S
s.S
K.dot(s.S)
K.dot(s.S)[:,0]
sum(K.dot(s.S)[:,0]*s.w)
sum(K.dot(s.S)[:,3]*s.w)
sum(K.dot(s.S)[:,4]*s.w)
sum(K.dot(s.S)[:,4])
K
phi = linalg.solve(K.dot(S), E[:,0])
phi = linalg.solve(K.dot(s.S), E[:,0])
phi
phi = linalg.solve(K, E[:,0])
linalg.norm(K.dot(S).dot(phi) - E[:,0])
linalg.norm(K.dot(s.S).dot(phi) - E[:,0])
phi = linalg.solve(K.dot(s.S), E[:,0])
linalg.norm(K.dot(s.S).dot(phi) - E[:,0])
K
phi
phi[0] = 0
linalg.norm(K.dot(s.S).dot(phi) - E[:,0])
K
max(abs(K))
abs(K)
help(map)
map(max, K)
K.dot(s.S).dot(phi)
linalg.norm(K.dot(s.S).dot(phi) - E[:,0])
phi = linalg.solve(K.dot(s.S), E[:,0])
linalg.norm(K.dot(s.S).dot(phi) - E[:,0])
K.dot(s.S).dot(phi)
phi
K.dot(s.S)
K.dot(s.S[:,1::])
linalg.norm(K.dot(s.S[:,1::]).dot(phi[1::]) - E[:,0])
linalg.norm(K.dot(s.S[:,0::]).dot(phi[0::]) - E[:,0])
K.dot(s.S[:,0])
s = sg.Segment(n, f_inargs = (sh.ellipse, (0, 3, 2)), quad='ps')
K
s.SL.T.dot(K)
reload(sg)
s = sg.Segment(n, f_inargs = (sh.ellipse, (0, 3, 2)), quad='ps')
s.SL.T.dot(K)
s.SL.T.dot(K)[0]
s.SL.T.dot(E[:,0])
phi1 = linalg.solve(K, E[:,0])
phi1 = linalg.solve(K.dot(s.S), E[:,0])
phi1
phi2 = linalg.solve(s.SL.T.dot(K.dot(s.S)), s.SL.T.dot(E[:,0]))
phi2
phi2 - phi1
phi3 = linalg.solve(s.SL.T.dot(K.dot(s.S))[1::, 1::], s.SL.T.dot(E[:,0])[1::])
phi3
phi3 -phi2[1::]
phi4 = linalg.solve(K.dot(s.S)[1::, 1::], E[1::,0])
phi4
