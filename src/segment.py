import matplotlib.pyplot as plt
import numpy as np
import cmath

pi = cmath.pi
cos = cmath.cos
sin = cmath.sin

def cZ(t, c = 0, a = 1, b = 1):
  return complex(c + a * cos(2 * pi * t) + 1j * b * sin(2 * pi * t))
def cZp(t, c = 0, a = 1, b = 1):
  return complex(2 * pi * a * cos(2 * pi * t) + 1j * 2 * pi * b *sin(2 * pi * t))
def cZpp(t, c = 0, a = 1, b = 1):
  return complex(2 * pi * 2 * pi * a * cos(2 * pi * t) + 1j * 2 * pi * 2 * pi * b * sin(2 * pi * t))

def sZ(t, p = 0, q = 1):
  return complex(p + t * (q - p))
def sZp(t, p = 0, q = 1):
  return complex(q - p)
def sZpp(t, p = 0, q = 1):
  return complex(0)

class Segment:
  def __init__(self, n, Z=cZ, Zp=cZp, Zpp=cZpp, args=[]):
    # trapezoidal nodes
    self.n = n
    self.t = np.array([float(k) / (n - 1) for k in range(n)], float)
    self.x = np.array([Z(t,*args) for t in self.t], np.cfloat)
    self.nx = np.array([-1j * Zp(t,*args) / abs(Zp(t,*args)) for t in self.t], np.cfloat)

def Z(t):
   return complex(cos(2 * pi * t) + 1j * sin(2 * pi * t))
def Zp(t):
  return complex(2 * pi * cos(2 * pi * t) + 1j * 2 * pi * sin(2 * pi * t))
def Zpp(t):
  return complex(2 * pi * 2 * pi * cos(2 * pi * t) + 1j * 2 * pi * 2 * pi * sin(2 * pi * t))

def circle(c = 0, r = 1):
  args=(c, r, r)
  return (cZ, cZp, cZpp, args)
def ellipse(c = 0, a = 1, b = 2):
  args=(c, a, b)
  return (cZ, cZp, cZpp, args)
def line(p = 0, q = 1):
  args=(p, q)
  return (sZ, sZp, sZpp, args)

(Z1, Zp1, Zpp1, args) = ellipse(1j, 5, 2)
(Z1, Zp1, Zpp1, args) = line(1j, 5)

s = Segment(n=10, Z=Z1, Zp=Zp1, Zpp=Zpp1, args=args)
plt.plot([x.real for x in s.x],[x.imag for x in s.x],'o-')
plt.axis('equal')
plt.show()


