import matplotlib.pyplot as plt
# import pylab
import numpy as np

def meshgrid((x1, x2, xn), (y1, y2, yn)=((), (), ())):
  xs = float(x2 - x1) / (xn - 1)
  x = np.arange(x1, x2 + xs, xs)
  if yn == ():
    (y1, y2, yn) = (x1, x2, xn)
  ys = float(y2 - y1) / (yn - 1)
  y = np.arange(y1, y2 + ys, ys)
  if x2 != x[-1]:
    print 'Warning: and of x meshgrid'
  if y2 != y[-1]:
    print 'Warning: and of y meshgrid'
  xx, yy = np.meshgrid(x, y, sparse=True)
  pp = xx + 1j *yy
  pp = pp.reshape(xn * yn)
  return (x, y, pp)

# z = np.sin(xx**2 + yy**2) / (xx**2 + yy**2)
# z = xx**2 + yy**2

def plot(h, x, y, vv, t='im'):
  v = vv.reshape((len(y), len(x)))
  if t == 'cf':
    h = plt.contourf(x, y, v, 20)
  elif t == 'im':
    h = plt.imshow(np.array(v, 'float64'), extent = (x[0], x[-1], y[0], y[-1]))
  # h = pylab.imshow(v, interpolation='nearest')
  plt.colorbar(h)
  plt.axis('equal')
  plt.show(block=False)
  return h
  
