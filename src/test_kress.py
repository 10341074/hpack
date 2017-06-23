import numpy as np

import layerpot as ly
import shapes as sh
import segment as sg
import mainpb as m

nso = 400
a = 3
b = 2

def test_kress_k(t, tau, a, b):
  return - a * b / (2 * np.pi) /( a**2 + b**2 - (a**2 - b**2) * np.cos(t + tau) )


so = sg.Segment(nso, f_inargs = (sh.ellipse, (0, a, b)), quad='p')
K = ly.layerpotSD_slf(s=so, t=sg.Pointset(x=[so.x[0]], nx=so.nx[0]), slf=1) 
i = sum(K[0])
# k = test_kress_k(so.t[0], so.t, a, b)

err = []
for n in range(4, 50):
  so = sg.Segment(n, f_inargs = (sh.ellipse, (0, a, b)), quad='ps')
  K = ly.layerpotSD_slf(s=so, t=sg.Pointset(x=[so.x[0]], nx=so.nx[0]), slf=1)
  err.append(sum(K[0]) - i)

