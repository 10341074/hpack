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
  plt.plot(np.log(x), np.log(y), 'k+-')
  plt.axis('equal')
  plt.axis('square')
  plt.title('Ciao')
  plt.show(block=False)
  return

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
  l = fig.levels
  # plt.colorbar()
  # p.plot_domain()
  path = fig.collections[1].get_paths()[0]
  v = path.vertices
  plt.show(block=False)
  return v

def stats(v):
  vc = np.array(v[:,0] + 1j * v[:, 1])
  bary = sum(vc) / v.size
  d = list(abs(vc - bary))
  ax1 = 2 * max(d)
  ax1ind = d.index(max(d))
  ax2 = 2 * min(d)
  ax2ind = d.index(min(d))
  return (bary, ax1, ax2, ax1ind, ax2ind)
