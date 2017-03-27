import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from matplotlib import cm

# import pylab
import numpy as np

# def meshgrid((x1, x2, xn), (y1, y2, yn)=((), (), ())):
def meshgrid(x1_x2_xn, y1_y2_yn=((), (), ())):
  x1, x2, xn = x1_x2_xn
  y1, y2, yn = y1_y2_yn
  xs = float(x2 - x1) / (xn - 1)
  x = np.concatenate( (np.arange(x1, x2, xs), [x2]) )
  if yn == ():
    (y1, y2, yn) = (x1, x2, xn)
  ys = float(y2 - y1) / (yn - 1)
  y = np.concatenate( (np.arange(y1, y2, ys), [y2]) )
  if x2 != x[-1]:
    print('Warning: end of x meshgrid')
  if y2 != y[-1]:
    print('Warning: end of y meshgrid')
  xx, yy = np.meshgrid(x, y, sparse=True)
  pp = xx + 1j *yy
  pp = pp.reshape(xn * yn)
  return (x, y, pp)

# z = np.sin(xx**2 + yy**2) / (xx**2 + yy**2)
# z = xx**2 + yy**2

def plot(x, y, vv, t='im', fig=[]):
  h = fig
  v = vv.reshape((len(y), len(x)))
  if t == 'cf':
    h = plt.contourf(x, y, v, 20)
    plt.colorbar()
  elif t == 'im':
    h = plt.imshow(np.array(v, 'float64'), extent = (x[0], x[-1], y[0], y[-1]))
    plt.colorbar()
  elif t == 'srf':
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    xx, yy = np.meshgrid(x, y, sparse=True)
    surf = ax.plot_surface(xx, yy, v, cmap=cm.coolwarm)
    fig.colorbar(surf)
    # surf = ax.plot_surface(X, Y, Z, cmap=cm.coolwarm, linewidth=0, antialiased=False)

  # h = pylab.imshow(v, interpolation='nearest')
  plt.axis('equal')
  plt.show(block=False)
  return h
  
