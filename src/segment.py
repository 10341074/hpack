#!../../venv/bin/python
import matplotlib.pyplot as plt
import numpy as np

import my_type
import plot
import layerpot as ly

from my_type import float_t
float = float_t
# print np.dtype(float)

# from mpl_toolkits.mplot3d import Axes3D
# from matplotlib import cm

pi = np.pi
cos = np.cos
sin = np.sin

def cZ(t, c = 0, a = 1, b = 1):
  return complex(c + a * cos(2 * pi * t) + 1j * b * sin(2 * pi * t))
def cZp(t, c = 0, a = 1, b = 1):
  return complex(-2 * pi * a * sin(2 * pi * t) + 1j * 2 * pi * b * cos(2 * pi * t))
def cZpp(t, c = 0, a = 1, b = 1):
  return complex(-2 * pi * 2 * pi * a * cos(2 * pi * t) - 1j * 2 * pi * 2 * pi * b * sin(2 * pi * t))

def sZ(t, p = 0, q = 1):
  return complex(p + t * (q - p))
def sZp(t, p = 0, q = 1):
  return complex(q - p)
def sZpp(t, p = 0, q = 1):
  return complex(0)

def kZ(t, a = []):
  return complex(cos(2*pi*t) + 0.65 * cos(2*2*pi*t) - 0.65 + 1j * 1.5 * sin(2*pi*t))
def kZp(t, a = []):
  return complex(-2*pi*sin(2*pi*t) - 2 * 0.65 * 2*pi*sin(2*2*pi*t) + 1j * 1.5 * 2*pi*cos(2*pi*t))
def kZpp(t, a = []):
  return complex(-2*pi*2*pi*cos(2*pi*t) - 2 * 2 * 0.65 * 2*pi*2*pi*cos(2*2*pi*t) - 1j * 1.5 * sin(2*pi*t))

class Segment:
  def __init__(self, n, periodic=False, Z=cZ, Zp=cZp, Zpp=cZpp, args=[], f=[], inargs=[]):
    if f != []:
      # print 'short form'
      (Z, Zp, Zpp, args) = f(*inargs)
    if periodic:
      self.t = np.array([float(k) / n for k in range(n)], float)
    else:
      self.t = np.array([float(k) / (n - 1) for k in range(n)], float)
    # trapezoidal nodes
    self.n = n
    self.x = np.array([Z(t, *args) for t in self.t], np.cfloat)

    self.dx = np.array([Zp(t, *args) for t in self.t], np.cfloat)
    self.speed = np.array(abs(self.dx), float)
    self.nx = np.array(-1j * self.dx / self.speed, np.cfloat)

    self.ddx = np.array([Zpp(t, *args) for t in self.t], np.cfloat)    
    # s.kappa = -real(conj(-1i*dZdt).*s.Zpp(s.t)) ./ s.speed.^3; %curvature
    self.kappa = -np.real(np.conj(-1j * self.dx) * self.ddx) / (self.speed**3) # signed curvature
    if periodic:
      self.w = np.array([sp / n for sp in self.speed], float)
    else:
      self.w = np.array([sp / (n-1) for sp in self.speed], float)
      self.w[0] = self.w[0] * 0.5
      self.w[-1] = self.w[-1] * 0.5
  def plot(self, p=False):
    xx = [x.real for x in self.x]
    yy = [x.imag for x in self.x]
    if p:
      xx.append(xx[0])
      yy.append(yy[0])
    plt.plot(xx, yy,'o-')


class Boundary:
  def __init__(self, pieces=[]):
    self.pc = pieces
    self.x = np.array([z for p in pieces for z in p.x])
    self.nx = np.array([nx for p in pieces for nx in p.nx])
    self.speed = np.array([sp for p in pieces for sp in p.speed])
    self.kappa = np.array([k for p in pieces for k in p.kappa])
    self.w = np.array([w for p in pieces for w in p.w])
  

class Pointset:
  def __init__(self, x=[]):
    self.x = np.array(x)

class Layer:
  def __init__(self, b=[], exp=[], dns=[]):
    self.b = b
    if dns == []:
      self.dns = np.ones(len(b.x))
    else:
      self.dns = dns
    self.exp = exp

def eval_layer(l, p, exp=[]):
  if exp == []:
    exp = l.exp
  A = exp(0, l.b, p, [])
  v = A.dot(l.dns)
  return v

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
###########################################################
if __name__ == "__main__":

  (Z1, Zp1, Zpp1, args1) = ellipse(1j, 5, 2)
  (Z2, Zp2, Zpp2, args2) = line(1j, 5)

  s = Segment(n=100, Z=Z1, Zp=Zp1, Zpp=Zpp1, args=args1)
  #plt.plot([x.real for x in s.x],[x.imag for x in s.x],'o-')
  #s = Segment(n=10, Z=Z2, Zp=Zp2, Zpp=Zpp2, args=args)
  #plt.plot([x.real for x in s.x],[x.imag for x in s.x],'o-')
  s = Segment(n=100, f=circle, inargs=(0,0.5), periodic=True)
  s = Segment(n=100, Z=kZ, Zp=kZp, Zpp=kZpp, args=[], periodic=True)

  z0 = 1j
  z1 = -2
  z2 = 2-0.5j
  z3 = 0.2j

  s1 = Segment(n=100, f=line, inargs=(z0, z1))
  s2 = Segment(n=100, f=line, inargs=(z1, z2))
  s3 = Segment(n=100, f=line, inargs=(z2, z3))
  s4 = Segment(n=100, f=line, inargs=(z3, z0))

  b = Boundary([s1, s2, s3, s4])
  b = Boundary([s])
  for s in b.pc:
    plt.plot([x.real for x in s.x],[x.imag for x in s.x],'o-')
    
  plt.axis('equal')
  plt.show()


  x, y, pp = plot.meshgrid((-2, 2, 100))
  pp = Pointset(pp)
  l = Layer(b, ly.layerpotS)


  vv = eval_layer(l, pp)
  # vv = vv.reshape((100,100))
  # h = plt.contourf(x,y,vv)
  # plt.colorbar(h)
  # plt.axis('equal')
  # plt.show()
  
  h = plt.figure()
  h = plot.plot(h, x,y,vv)
  
  plt.plot([x.real for x in s.x],[x.imag for x in s.x],'o-')
  plt.show()
  
  
  # fig = plt.figure()
  # ax = fig.add_subplot(111, projection='3d')
  
  # ax.plot_surface(xx,yy,vv)
  # surf = ax.plot_surface(xx, yy,vv, rstride=1, cstride=1, cmap=cm.coolwarm, linewidth=0, antialiased=False)
  # fig.colorbar(surf, shrink=0.5, aspect=5)
  # ax.set_xlabel('X Label')
  # ax.set_ylabel('Y Label')
  # ax.set_zlabel('Z Label')
  
  # plt.show()


  # plt.figure(1)
  # plt.imshow(your_first_image)
  
  # plt.figure(2)
  # plt.imshow(your_second_image)

  # plt.show(block=False) # That's important 

  # raw_input("Press ENTER to exist")
  
  vv = eval_layer(l, [])
  h = plt.figure()
  
  plt.plot(range(len(vv)),np.real(vv),'o-')
  plt.show()
