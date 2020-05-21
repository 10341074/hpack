import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from matplotlib import cm

import numpy as np

def meshgrid(x1_x2_xn, y1_y2_yn=((), (), ())):
  '''
  This function requires a triplet (x1, x2, xn), it's used for y if not specified
  It returns (x, y, pp), where 'pp' is the dyadic product of the axes, it's complex and flatten for easy evaluation
  The output (x, y, pp) is ready for plot(x, y, vv=fun(pp)) evaluated on 'pp'
  '''
  x1, x2, xn = x1_x2_xn
  y1, y2, yn = y1_y2_yn
  xs = float(x2 - x1) / (xn - 1)
  x = np.concatenate( ([x1 + k * xs for k in range(xn-1)], [x2]) )
  if yn == ():
    (y1, y2, yn) = (x1, x2, xn)
  ys = float(y2 - y1) / (yn - 1)
  y = np.concatenate( ([y1 + k * ys for k in range(yn-1)], [y2]) )
  if x2 != x[-1]:
    print('Warning: end of x meshgrid')
  if y2 != y[-1]:
    print('Warning: end of y meshgrid')
  xx, yy = np.meshgrid(x, y, sparse=True)
  pp = xx + 1j *yy
  pp = pp.reshape(xn * yn)
  return (x, y, pp)

def plot(x, y, vv, show=1, pltype='im', pllevels='exp', colorbar=True, *args, **kwargs):
  fig = plt.figure()
  v = vv.reshape((len(y), len(x)))
  if pltype == 'c':
    if 'levels' not in kwargs: kwargs['levels'] = 30
    if pllevels == 'log': kwargs['levels'] = (np.log(np.linspace(1, 2, kwargs['levels'])) / np.log(2)) * (np.max(v) - np.min(v)) + np.min(v)
    if pllevels == 'exp': kwargs['levels'] = (np.exp(np.linspace(0, 1, kwargs['levels'])) - 1) / (np.exp(1) - 1) * (np.max(v) - np.min(v)) + np.min(v)
    fig = plt.contour(x, y, v, *args, **kwargs)
    if colorbar: plt.colorbar()
  if pltype == 'cf':
    if 'levels' not in kwargs: kwargs['levels'] = 30
    if pllevels == 'log': kwargs['levels'] = (np.log(np.linspace(1, 2, kwargs['levels'])) / np.log(2)) * (np.max(v) - np.min(v)) + np.min(v)
    if pllevels == 'exp': kwargs['levels'] = (np.exp(np.linspace(0, 1, kwargs['levels'])) - 1) / (np.exp(1) - 1) * (np.max(v) - np.min(v)) + np.min(v)
    fig = plt.contourf(x, y, v, *args, **kwargs)
    if colorbar: plt.colorbar()
  elif pltype == 'im':
    fig = plt.imshow(np.array(v[::-1], 'float64'), extent = (x[0], x[-1], y[0], y[-1]), *args, **kwargs)
    if colorbar: plt.colorbar()
  elif pltype == 'srf' or pltype == 'maxsrf':
    if pltype == 'maxsrf':
      mng = plt.get_current_fig_manager()
      mng.resize(*mng.window.maxsize())
    ax = fig.gca(projection='3d')
    xx, yy = np.meshgrid(x, y,  sparse=True)
    surf = ax.plot_surface(xx, yy, v, cmap=cm.coolwarm)
    if colorbar: fig.colorbar(surf)
  # ------------------
    # surf = ax.plot_surface(X, Y, Z, cmap=cm.coolwarm, linewidth=0, antialiased=False)
    # h = pylab.imshow(v, interpolation='nearest')
  # ----------------------------------------------------------------------------------------------------
  plt.axis('square')
  plt.axis('equal')
  if show:
    plt.show(block=False)
  return fig

def contour(x, y, vv, val=1e-3):
  fig = plt.figure()
  v = vv.reshape((len(y), len(x)))
  fig = plt.contour(x, y, v, [val])
  return fig.collections[0].get_paths()[0].vertices

def surf(xx, yy, vv):
  fig = plt.figure()
  ax = fig.gca(projection='3d')
  # xx, yy = np.meshgrid(x, y, sparse=True)
  surf = ax.plot_surface(xx, yy, vv, cmap=cm.coolwarm)
  fig.colorbar(surf)
  plt.axis('equal')
  plt.show(block=False)
  return fig
