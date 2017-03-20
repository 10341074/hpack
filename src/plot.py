import matplotlib.pyplot as plt
import numpy as np

def meshgrid((x1, x2, xn), (y1, y2, yn)=((), (), ())):
  xs = float(x2 - x1) / xn
  x = np.arange(x1, x2, xs)
  if yn == ():
    (y1, y2, yn) = (x1, x2, xn)
  ys = float(y2 - y1) / yn
  y = np.arange(y1, y2, ys)
  xx, yy = np.meshgrid(x, y, sparse=True)
  pp = xx + 1j *yy
  pp = pp.reshape(xn * yn)
  return (x, y, pp)

# z = np.sin(xx**2 + yy**2) / (xx**2 + yy**2)
# z = xx**2 + yy**2

def plot(h, x, y, vv):
  v = vv.reshape((len(y), len(x)))
  np.shape(v)
  h = plt.contourf(x, y, v, 20)
  plt.colorbar(h)
  plt.axis('equal')
  plt.show(block=False)
  return h
