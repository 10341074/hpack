# coding: utf-8
import my_session_150417_1200
s
from my_session_150417_1200 import *
s
Sb = linf.gramschmidtBoundary(s)
Sb
Sb.T.dot(Sb)
Sb.T.dot(np.diagflat(s.w).dot(Sb))
S
Sb
L0 = ipb.computeallpsiL0(s, E)
L0b = ipb.computeallpsiL0(s, Sb[:, 1::])
e
eb = linalg.eigvals(L0b)
plt.plot(e.real, e.imag, 'bo')
plt.show(block=False)
e = linalg.eigvals(L0)
plt.plot(e.real, e.imag, 'bo')
plt.show(block=False)
plt.plot(e2.real, e2.imag, 'ro')
plt.show(block=False)
eb
e2
e2 = linalg.eigvals(L02)
eb 
plt.plot(eb.real, eb.imag, 'go')
plt.show(block=False)
e2
max(e2)
S.T.dot(S)
L02 = ipb.computeallpsiL0(s, S)
e2 = linalg.eigvals(L02)
max(e2)
e2
max(eb)
eb = linalg.eigvals(L0b)
L0b.shape
L02.shape
s.w
