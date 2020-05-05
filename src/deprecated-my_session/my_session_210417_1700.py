# coding: utf-8
from __mainpb__ import *
LL0 = ipb.computeLL0(ld, so, T=np.eye(so.n), c=c, testset=0)
LL0B = ipb.computeLL0B(ld, so, T=so.B[:, 1:], c=c, testset=0)
LL0.shape
LL0B.shape
(w, v) = linalg.eig(LL0B)
plt.plot(w.real, w.imag, 'b*')
plt.show()
LLdiff = ipb.computeLLdiff(ld, so, np.eye(so.n), c)
LLBdiff = ipb.computeLLBdiff(ld, so, np.eye(so.n), c)
LLdiff.shape
LLBdiff.shape
LLBdiff = ipb.computeLLBdiff(ld, so, so.B[:, 1:], c)
LLBdiff.shape
(w, v) = linalg.eig(LLBdiff)
plt.plot(w.real, w.imag, 'b*')
plt.show()
v
issym(LLdiff)
(w, v) = linalg.eig(0.5*(LLBdiff + LLBdiff.T))
plt.plot(w.real, w.imag, 'b*')
plt.show()
w
plt.plot(w.real, w.imag, 'b*')
plt.show()
(w, v) = linalg.eig(0.5*(LLdiff + LLdiff.T))
plt.plot(w.real, w.imag, 'b*')
plt.show()
T = np.eye(so.n)
allpsi0 = dpb.mapNtoD0(so, T, so.s0) 
Lo = ly.layerpotSD(s=so, t =ld)
rhsdiff = Lo.dot(allpsi0)
allpsi = dpb.mapNtoDdiff(so, ld, rhsdiff, c, so.s0)
x, y, p = ipb.meshgrid((-5, 5, 50))
x, y, p = ipb.meshgrid((-8, 8, 50))
x, y, p = ipb.meshgrid((-8, 8, 60))
T = v[:, 22]
allpsi0 = dpb.mapNtoD0(so, T, so.s0) 
Lo = ly.layerpotSD(s=so, t =ld)
rhsdiff = Lo.dot(allpsi0)
allpsi = dpb.mapNtoDdiff(so, ld, rhsdiff, c, so.s0)
u = ly.layerpotS(s=so, t=p).dot(allpsi[0:so.n]) + ly.layerpotS(s=ld, t=p).dot(allpsi[-ld.n::])
plot.plot(x,y,u)
allpsi
v[:,22]
T = v[:, 22].real
allpsi0 = dpb.mapNtoD0(so, T, so.s0) 
Lo = ly.layerpotSD(s=so, t =ld)
rhsdiff = Lo.dot(allpsi0)
allpsi = dpb.mapNtoDdiff(so, ld, rhsdiff, c, so.s0)
u = ly.layerpotS(s=so, t=p).dot(allpsi[0:so.n]) + ly.layerpotS(s=ld, t=p).dot(allpsi[-ld.n::])
plot.plot(x,y,u)
x, y, p = ipb.meshgrid((-2, 2, 60))
u = ly.layerpotS(s=so, t=p).dot(allpsi[0:so.n]) + ly.layerpotS(s=ld, t=p).dot(allpsi[-ld.n::])
plot.plot(x,y,u)
sum(v[:, 22].real**2*so.w)
v.shape
LLdiff.dot(v[:, 22].real / 0.23)
v[:, 22]
sd = sg.Segment(nsd, Z_args = (sh.kZ, sh.kZp, sh.kZpp, ()), quad='gp', aff=(0, 0.8 + 0.8j))
sd = sg.Segment(nsd, Z_args = (sh.kZ, sh.kZp, sh.kZpp, ()), quad='gp', aff=(0, 0.8 + 0.8j))
reload(sh)
sd = sg.Segment(nsd, Z_args = (sh.kZ, sh.kZp, sh.kZpp, ()), quad='gp', aff=(0, 0.8 + 0.8j))
sd.plot()
plt.show
plt.show()
sd = sg.Segment(nsd, Z_args = (sh.kZ, sh.kZp, sh.kZpp, ()), quad='ps', aff=(0, 0.8 + 0.8j))
sd.plot()
plt.show()
