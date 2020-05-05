# coding: utf-8
# check on eigenvalues
from __load__ import *
s = sg.Segment(40, f_inargs=(sh.ellipse, (0, 3, 2)), quad='ps')
s.plot()
plt.show()
import linfunc as linf
import linfunc as linf
reload(linf)
S = linf.gramschmidt(np.ones(40))
S.T.dot(S)
E = np.diagflat(np.ones(40))
E
linalg
help(linalg)
e1 =linalg.eig(E)
e1
e1 =linalg.eigvals(E)
e1
e1 =linalg.eigvals(S)
e1
max(abs(e1 - np.ones(40)))
abs(e1 - np.ones(40))
E.shape
S.shape
L = S.T.dot(E.dot(S))
L
e1 =linalg.eigvals(L)
e1
L0 = ipb.computeallpsiL0(s, E)
L02 = ipb.computeallpsiL0(s, S)
e = linalg.eigvals(L0)
e
plt.plot(e.real, e.imag, 'bo')
plt.show(block=False)
e2 = linalg.eigvals(L02)
plt.plot(e2.real, e2.imag, 'ro')
plt.show(block=False)
plt.plot(e.real, e.imag, 'bo')
plt.show(block=False)
plt.plot(e.real, e.imag, 'bo')
plt.show(block=False)
plt.figure()
plt.plot(e2.real, e2.imag, 'ro')
plt.show(block=False)
reload(ipb)
plt.figure()
L0 = ipb.computeallpsiL0(s, E)
L02 = ipb.computeallpsiL0(s, E)
e = linalg.eigvals(L0)
e2 = linalg.eigvals(L02)
plt.plot(e.real, e.imag, 'bo')
plt.show(block=False)
plt.show(block=False)
e2 = linalg.eigvals(L02)
plt.show(block=False)
plt.figure()
plt.show(block=False)
plt.show(block=False)
plt.figure()
plt.show(block=False)
plt.plot(e.real, e.imag, 'bo')
plt.show(block=False)
plt.figure()
plt.plot(e2.real, e2.imag, 'ro')
plt.show(block=False)
# ok same plots
plt.plot(e2.real, e2.imag, 'ro')
plt.show(block=False)
L03 = ipb.computeallpsiL0(s, S[:, 1::]) # not Boundary, only without one eig: same eigs?
e3 = linalg.eigvals(L03)
plt.plot(e3.real, e3.imag, 'ro')
plt.show(block=False)
plt.plot(e2.real, e2.imag, 'bo')
plt.show(block=False)
plt.plot(e3.real, e3.imag, 'ro')
plt.show(block=False)
# get_ipython().magic('save my-session-150417-1200 1-78')
plt.close("all")
