from src import *
import numpy as np
import scipy.linalg as linalg

import layerpot as ly
import shapes as sh
import segment as sg
import geometries as gm

import mainpb as m

import inverseproblem as ipb

def v_ex(z):
  return z.real**3 - 3 * z.real * z.imag**2
def v_p_ex(z):
  return (3 * z.real**2 - 3 * z.imag**2) + 1j * ( - 6 * z.real * z.imag)

def e0(tt, st0, st1):
  return (tt - st0) / (st1 - st0)
def e1(tt, st0, st1):
  return (st1 - tt) / (st1 - st0)
def truncate(f):
  return f * (f>=0) * (f<=1)
def truncate_up(f):
  return 1 + (f - 1) * (f<=1)
def elem_deg1(s, t):
  st_ext = np.concatenate( ([s.t[-1] - 1], s.t, [s.t[0] + 1]) )
  A = np.empty((len(t.t), len(s.t)))
  for j in range(1, len(st_ext) - 1 ):
    A[:, j-1] = truncate_up(truncate(e0(t.t, st_ext[j-1], st_ext[j])) + truncate(e1(t.t, st_ext[j], st_ext[j+1])))\
              + truncate_up(truncate(e0(t.t - 1, st_ext[j-1], st_ext[j])) + truncate(e1(t.t - 1, st_ext[j], st_ext[j+1])))\
              + truncate_up(truncate(e0(t.t + 1, st_ext[j-1], st_ext[j])) + truncate(e1(t.t + 1, st_ext[j], st_ext[j+1])))
  return A  

def iters(rng, f_ex, s_ex):
  err = []
  for n in rng:
    f, s = solve(n)
    # A = elem_deg1(s, s_ex)
    # new_err = linalg.norm((A.dot(f) - f_ex) * s_ex.w)
    f_ex_n = v_ex(s.x)
    f_ex_n = f_ex_n - sum(f_ex_n * s.w) / sum(s.w)
    new_err = linalg.norm((f - f_ex_n) * s.w)
    err = np.concatenate((err, [new_err]))
    # plt.plot(s_ex.t, f_ex, '+-')
    # plt.plot(s.t, f, '+-')
    # plt.plot(s_ex.t, A.dot(f),'+-')
    # plt.show(block=False)
    # end = input('Press')
  return err

def solve(n):
  s = gm.sg_one_kite(n)
  # s = gm.sg_one_triangle(n)
  # L0 = ipb.computeL0(so = s, T = np.eye(s.n))
  v_p = ly.scalar(v_p_ex(s.x), s.nx)
  v = ipb.computeL0(so = s, T = v_p)
  # return L0.dot(v_p)
  return v, s

def exact(n):
  s_ex = gm.sg_one_kite(n)
  f_ex = v_ex(s_ex.x)
  f_ex = f_ex - sum(f_ex * s_ex.w) / sum(s_ex.w)
  return f_ex, s_ex
  

def total(n_ex):
  f_ex, s_ex = exact(n_ex)
  err = iters(range(10, n_ex, 10), f_ex, s_ex)
  plt.plot(err, '+-')
  plt.show(block=False)
  return err
