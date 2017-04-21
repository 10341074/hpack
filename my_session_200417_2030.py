# coding: utf-8
from __mainpb__ import *
LL0B = ipb.computeLL0B(ld, so, sb, c, testset=0)
LB = ipb.computeLB(ld, so, sb, c)
LB = ipb.computeLB(so, ls, c)
LB = ipb.computeLB(so, ld, c)
LB.shape
L0B = ipb.computeL0B(so)
(w, v) = linalg.eig(LL0B)
min(w)
(w, v) = linalg.eig(LL0B)
w[22]
LL0B.dot(v[22])
f=LL0B.dot(v[22])
f-v[22]
help(linalg.eig)
f=LL0B.dot(v[:,22])
f-v[22]
f-v[:,22]
help(linalg.eig)
(w, v) = linalg.eig(LL0B)
w.shape
v.shape
f=LL0B.dot(v[:, 22])
w[22]
ff = w[22] * v[:, 22]
f -ff
w[22]
v[:, 22]
f
f.shape
f1 = L0B.dot(f)
f1 = L0B.dot(v[:, 22])
f2 = LB.dot(v[:, 22])
f1 - f2
f1
f2
L = ipb.computeL0(so, so.B[:, 1:])
L0 = ipb.computeL0(so, so.B[:, 1:])
L = ipb.computeL(so, ld, so.B[:, 1:], c)
L.shape
L0.shape
g = L.dot(v[:,22])
g0 = L0.dot(v[:, 22])
g- g0
f = v[:, 22]
psi0 =  dpb.mapNtoD0(so, f, so.s0)
so.plot()
plt.show()
so.s0
f= np.concatenate(( [0], f))
f = so.B[:, 1::].dot(v[:, 22])
vf = so.B[:, 1::].T.dot(np.diagflat(so.w).dot(f))
vf
vf = so.B.T.dot(np.diagflat(so.w).dot(f))
vf
vf[1:]-v[:,22]
max(abs(vf[1:]-v[:,22]))
f.shape
psi0 = dpb.mapNtoD0(so, f, so.s0)
psi = dpb.mapNtoD(so, ld, f, c)
psi.shape
ld.n
q = range(8)
q[0:6]
list(q[0:6])
d0 = ly.layerpotS(s=so).dot(psi)
d0 = ly.layerpotS(s=so).dot(psi0)
d = ly.layerpotS(s=so).dot(psi[0:so.n]) + ly.layerpotS(s=ld,t=so).dot(psi[-ld.n::])
d0
d
d-d0
sum(d * so.w)
x
x, y, p = ipb.meshgrid(( - 8, 8, 50))
x, y, p = ipb.meshgrid(( -8, 8, 50))
x, y, p = ipb.meshgrid(( -8, 8, 40))
x, y, p = ipb.meshgrid(( -8, 8, 60))
u = ly.layerpotS(s=so, t = p).dot(psi0)
plot(x, y, u)
plot.plot(x, y, u)
f
sum(f)
psi0
plt.plot(range(so.n), psi0.real, 'b*')
plt.show()
plt.plot(range(so.n), psi0.real, 'b*')
plt.show()
plt.plot(range(so.n), psi0.real, 'b*')
plt.plot(range(so.n), psi0.real, 'b-*')
plt.show()
psi0 = dpb.mapNtoD0(so, np.eye(so.n)[:,0], so.s0)
u = ly.layerpotS(s=so, t = p).dot(psi0)
plt.plot(range(so.n), psi0.real, 'b*')
plt.show()
psi0
plot.plot(x, y, u)
plot.plot(x, y, u,'srf')
reload(ipb)
reload(ipb)
reload(ipb)
LL0 = ipb.computeLL0(ld, so, np.eye(so.n), c)
LL0.shape
(w2, v2) = linalg.eig(LL0)
plt.plot(w2.real, w2.imag, 'b*')
plt.show()
w2.shape
plt.show()
plt.plot(w2.real, w2.imag, 'b*')
plt.show()
w2
psi0
psi0 = dpb.mapNtoD0(so, f, so.s0)
sum(f**2*so.w)
LL0B
issym(LL0B)
LL0Bs = (LL0B + LL0B.T) / 0.5
(ww, vv) = linalg.eig(LL0Bs)
ww
min(ww)
max(ww)
d0
d
plt.plot(range(so.n), psi.real, 'b*')
plt.plot(range(so.n+ ld.n), psi.real, 'b*')
plt.show()
plt.plot(range(ld.n), psi[-ld.n::].real, 'b*')
plt.show()
u1 = ly.layerpotS(s=so,t=p).dot(psi[0:so.n]) + ly.layerpotS(s=ld,t=p).dot(psi[-ld.n::])
plot.plot(x, y, u1, 'srf')
u1.shape
x
u1
psi
psi
psi0
f
psi0 = dpb.mapNtoD0(so, f.real, so.s0)
psi = dpb.mapNtoD(so, ld, f.real, c)
psi = dpb.mapNtoD(so, ld, f.real, c, so.s0)
psi
u1 = ly.layerpotS(s=so,t=p).dot(psi[0:so.n]) + ly.layerpotS(s=ld,t=p).dot(psi[-ld.n::])
plot.plot(x, y, u1, 'srf')
plot.plot(x, y, u1, 'srf')
d = ly.layerpotS(s=so).dot(psi[0:so.n]) + ly.layerpotS(s=ld,t=so).dot(psi[-ld.n::])
d0 = ly.layerpotS(s=so).dot(psi0)
d-d0
psi0
plt.plot(range(so.n), psi0, 'b*')
plt.show()
d0
plt.plot(range(so.n), d0, 'b*')
plt.show()
plt.plot(range(so.n+ld.n), psi, 'b*')
plt.show()
plt.plot(range(so.n+ld.n), d, 'b*')
plt.plot(range(so.n), d, 'b*')
plt.show()
plt.plot(range(so.n), d, 'b-*')
plt.plot(range(so.n), d0, 'r-*')
plt.show()
f
K=ly.layerpotSD(s=so)+0.5*np.eye(so.n)
K.dot(psi0) - f
max(abs(K.dot(psi0) - f))
Ko=ly.layerpotSD(s=so)+0.5*np.eye(so.n)
Ko.dot(psi[0:so.n]) + ly.layerpotSD(s=ld, t=so).dot(psi[-ld.n::]) - f
max(abs(Ko.dot(psi[0:so.n]) + ly.layerpotSD(s=ld, t=so).dot(psi[-ld.n::]) - f))
ly.layerpotSD(s=so, t=ld).dot(psi[0:so.n]) + ly.layerpotSD(s=ld).dot(psi[-ld.n::]) + c*0.5*psi[-ld.n:]
c
psi
ly.layerpotSD(s=so, t=ld).dot(psi[0:so.n])
ww
plt.plot(range(so.n), w, 'r-*')
plt.plot(range(len(w)), w, 'r-*')
plt.show()
LL0B = ipb.computeLL0B(ld, so, (), c)
sd = sg.Segment(nsd, f_inargs = (sh.ellipse, (0, rd, 1.5*rd)), quad='ps')
sd = sg.Segment(40, f_inargs = (sh.ellipse, (0, rd, 1.5*rd)), quad='ps')
sd = sg.Segment(40, f_inargs = (sh.ellipse, (0, 1, 1.5)), quad='ps')
so = sg.Segment(140, f_inargs = (sh.circle, (0, 6)), quad='ps')
LL0B = ipb.computeLL0B(sd, so, (), c)
LL0Bs = (LL0B + LL0B.T) / 0.5
(ww, vv) = linalg.eig(LL0Bs)
plt.plot(range(len(w)), w, 'r-*')
plt.show()
sd = sg.Segment(80, f_inargs = (sh.ellipse, (0, 1, 1.5)), quad='ps')
sd.plot()
so.plot()
plt.show()
LL0B = ipb.computeLL0B(sd, so, (), c)
LL0Bs = (LL0B + LL0B.T) / 0.5
(ww, vv) = linalg.eig(LL0Bs)
plt.plot(range(len(w)), w, 'r-*')
plt.show()
so = sg.Segment(140, f_inargs = (sh.circle, (0, 3)), quad='ps')
so.plot()
sd.plot()
plt.show()
so = sg.Segment(180, f_inargs = (sh.circle, (0, 2.5)), quad='ps')
so.plot()
sd.plot()
plt.show()
LL0B = ipb.computeLL0B(sd, so, (), c)
LL0B.shape
LL0Bs = (LL0B + LL0B.T) / 0.5
(ww, vv) = linalg.eig(LL0Bs)
plt.plot(range(len(w)), w, 'r-*')
plt.show()
len(ww)
plt.plot(range(len(w)), w, 'r-*')
plt.show()
plt.plot(range(len(ww)), ww, 'r-*')
plt.show()
so = sg.Segment(80, f_inargs = (sh.circle, (0, 2.5)), quad='ps')
sd = sg.Segment(80, f_inargs = (sh.ellipse, (0, 1, 1.5)), quad='ps')
LL0B = ipb.computeLL0B(sd, so, (), c)
LL0Bs = (LL0B + LL0B.T) / 0.5
(ww, vv) = linalg.eig(LL0Bs)
plt.plot(range(len(ww)), ww, 'r-*')
plt.show()
sd = sg.Segment(120, f_inargs = (sh.ellipse, (0, 1, 1.5)), quad='ps')
so = sg.Segment(120, f_inargs = (sh.circle, (0, 2.5)), quad='ps')
LL0B = ipb.computeLL0B(sd, so, (), c)
LL0Bs = (LL0B + LL0B.T) / 0.5
(ww, vv) = linalg.eig(LL0Bs)
plt.plot(range(len(ww)), ww, 'r-*')
plt.show()
vv.T.dot(vv)
plt.plot(range(len(ww)), np.log(ww), 'r-*')
plt.show()
ww
ww.real
plt.plot(range(len(ww)), np.log(ww.real), 'r-*')
plt.show()
plt.plot(range(len(ww)), np.log(ww.real), 'r-*')
plt.plot(range(len(ww)), np.log(ww), 'r-*')
plt.show()
np.log(ww)
vv
plt.plot(range(len(ww)), np.log(ww.real[0:20]), 'r-*')
plt.plot(range(len(m)), np.log(ww.real[0:m]), 'r-*')
m = 20
plt.plot(range(len(m)), np.log(ww.real[0:m]), 'r-*')
plt.plot(range(m), np.log(ww.real[0:m]), 'r-*')
m = 10
plt.plot(range(m), np.log(ww.real[0:m]), 'r-*')
m = 14
m = 6
plt.plot(range(m), np.log(ww.real[0:m]), 'r-*')
m = 4
plt.plot(range(m), np.log(ww.real[0:m]), 'r-*')
m = 3
plt.plot(range(m), np.log(ww.real[0:m]), 'r-*')
ww.real[0:3]
plt.plot(range(m), np.log(-ww.real[0:m]), 'r-*')
m = 20
plt.plot(range(m), np.log(-ww.real[0:m]), 'r-*')
m = 40
plt.plot(range(m), np.log(-ww.real[0:m]), 'r-*')
ww.shape
m = 80
plt.plot(range(m), np.log(-ww.real[0:m]), 'r-*')
m = 60
plt.plot(range(m), np.log(-ww.real[0:m]), 'r-*')
m = 50
plt.plot(range(m), np.log(-ww.real[0:m]), 'r-*')
m = 40
plt.plot(range(m), np.log(-ww.real[0:m]), 'r-*')
plt.show()
plt.plot(range(m), np.log(-ww.real[0:m]), 'r-*')
plt.show()
-ww.real[0:50]
for k in range(len(ww)):
    if -ww.real < 0:
        print(k)
        
for k in range(len(ww)):
    if -ww[k].real < 0:
        print(k)
        
