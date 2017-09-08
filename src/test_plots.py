from src import *
import numpy as np
import scipy.linalg as linalg

import layerpot as ly
import shapes as sh
import segment as sg
import geometries as gm

import mainpb as m

import inverseproblem as ipb

def plot_loglog(x=(), y=()):
  if x==(): x=range(1, len(y) + 1)
  fig = plt.figure()
  plt.plot(np.log(x), np.log(y), 'k+-', lw=1, ms=4)
  plt.axis('equal')
  plt.axis('square')
  plt.show(block=False)
  return fig
def plot_loglogscale(x=(), y=()):
  fig = plt.figure()
  plt.plot(x, y, 'k+-', lw=1, ms=4, ls=':')
  ax = fig.add_subplot(111)
  ax.set_yscale('log')
  ax.set_xscale('log')
  plt.show(block=False)
  return fig

def grid(x, y):
  for xk in x:
    plt.plot([xk, xk], [y[0], y[-1]], 'k', linestyle='--')
  for yk in y:
    plt.plot([x[0], x[-1]], [yk, yk], 'k', linestyle='--')
  return

def plot_contourf(p, args):
  z = np.array(p.z)
  for k in range(p.x.size):
    if p.p.flag_inside_s[k] == 0:
      z[k] = 0
  z = z.reshape((p.x.size, p.y.size))
  ##################
  fig = plt.figure()
  fig = plt.contourf(p.x, p.y, z, 20)
  l = fig.levels
  # plt.colorbar()
  # p.plot_domain()
  plt.show(block=False)
  # plt.clf()
  ##################
  fig = plt.figure()
  l1 = l[0] + 0.05*(l[-1] - l[0])
  fig = plt.contourf(p.x, p.y, z, [l[0], l1 , l[-1]])
  plt.colorbar()
  p.plot_domain()
  
  plt.show(block=False)
  ##################
  # fig = plt.figure()
  # z = z.reshape(p.x.size*p.y.size)
  # for k in range(p.x.size * p.y.size):
  #   if z[k] > l1:
  #     z[k] = 0
  # z = z.reshape((p.x.size, p.y.size))
  # fig = plt.contourf(p.x, p.y, z, 20)
  # plt.colorbar()
  # p.plot_domain()
  # plt.show(block=False)
  #################
  # fig = plt.figure()
  # z = np.array(p.z)
  # for k in range(p.x.size * p.y.size):
  #   if z[k] < l1:
  #     z[k] = l1
  # z =p.z.reshape((p.x.size, p.y.size))
  # fig = plt.contourf(p.x, p.y, z, 20, vmin = l1, vmax=l[-1])
  # plt.colorbar()
  # p.plot_domain()
  # plt.show(block=False)
  return
  
# p = fig.collections[0].get_paths()[0]
# v = p.vertices
# x = v[:,0]
# y = v[:,1]
def plot_contourf_1(p, args):
  z = np.array(p.z)
  for k in range(p.x.size):
    if p.p.flag_inside_s[k] == 0:
      z[k] = 0
  z = z.reshape((p.x.size, p.y.size))
  fig = plt.figure()
  fig = plt.contourf(p.x, p.y, z, 1)
  plt.axis('equal')
  plt.axis('square')
  # ax = fig.add_subplot(111)
  # ax.set_xlim(-3, 3)
  # ax.set_ylim(-3, 3)
  plt.show(block=False)
  l = fig.levels
  # plt.colorbar()
  # p.plot_domain()
  # path = fig.collections[1].get_paths()[1] # triangleplot
  path = fig.collections[1].get_paths()[0]
  v = path.vertices
  # v = np.concatenate(( v[-1:-1], v[:-6])) # triangleplot
  plt.show(block=False)
  return v

def stats(v):
  vc = np.array(v[:, 0] + 1j * v[:, 1])
  # bary = sum(vc) / v.size
  vtot = np.concatenate((v[:,0], v[:,1]))
  x1_x2 = ( min(np.floor(vtot)), max(np.ceil(vtot)) )
  bary = barycenter(x1_x2, 5000, v[:, 0], v[:, 1])
  d = list(abs(vc - bary))
  ax1 = max(d)
  ax1ind = d.index(max(d))
  ax2 = min(d)
  ax2ind = d.index(min(d))
  return (bary, ax1, ax2, ax1ind, ax2ind)

# https://stackoverflow.com/questions/217578/how-can-i-determine-whether-a-2d-point-is-within-a-polygon
def pnpoly(vertx, verty, testx, testy):
  nvert = vertx.size
  i, j, c = 0, nvert-1, 0
  for i in range(nvert):
    if ( ((verty[i]>testy) != (verty[j]>testy)) and (testx < (vertx[j]-vertx[i]) * (testy-verty[i]) / (verty[j]-verty[i]) + vertx[i]) ):
      c = 1 - c
    j = i
    i = i + 1
  return c

def barycenter(x1_x2, xn, vertx, verty):
  x1, x2 = x1_x2
  testx = numpy.random.uniform(x1, x2, xn)
  testy = numpy.random.uniform(x1, x2, xn)
  bary, count = 0, 0
  for i in range(xn):
    if pnpoly(vertx, verty, testx[i], testy[i]) == 1:
      bary = bary + testx[i] + 1j * testy[i]
      count = count + 1
  return bary / count
