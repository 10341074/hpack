import numpy as np
import matplotlib.pyplot as plt
import shapes as sh
import segment as sg

###########################################################################################################
def e0(tt, st0, st1):
  # first element function
  return (tt - st0) / (st1 - st0)
def e1(tt, st0, st1):
  # second element function
  return (st1 - tt) / (st1 - st0)
def truncate(f):
  return f * (f>=0) * (f<=1)
def truncate_up(f):
  # truncate in the exact quadrature nodes the sum of two element functions
  return 1 + (f - 1) * (f<=1)
def elem_deg1(s, t):
  # compute all element functions for sources s and targets t
  st_ext = np.concatenate( ([s.t[-1] - 1], s.t, [s.t[0] + 1]) )
  A = np.empty((len(t.t), len(s.t)))
  for j in range(1, len(st_ext) - 1 ):
    A[:, j-1] = truncate_up(truncate(e0(t.t, st_ext[j-1], st_ext[j])) + truncate(e1(t.t, st_ext[j], st_ext[j+1])))\
              + truncate_up(truncate(e0(t.t - 1, st_ext[j-1], st_ext[j])) + truncate(e1(t.t - 1, st_ext[j], st_ext[j+1])))\
              + truncate_up(truncate(e0(t.t + 1, st_ext[j-1], st_ext[j])) + truncate(e1(t.t + 1, st_ext[j], st_ext[j+1])))
  return A  
###########################################################################################################

def normErr(sFun, nbig, nsml, g, gh, sbig=(), ssml=(), norm = 'L2'):
  if sbig == ():
    sbig = sFun(nbig)
  if ssml == ():
    ssml = sFun(nsml)
  A = elem_deg1(ssml, sbig)
  g_ex = g(sbig.x); g_ex = g_ex - sum(g_ex * sbig.w) / sum(sbig.w)
  g_n = A.dot(gh); g_n = g_n - sum(g_n * sbig.w) / sum(sbig.w)
  if norm == 'L2':# metti media zero??
    err = np.sqrt(sum(sbig.w * (g_ex - g_n)**2))
  elif norm == 'Inf':
    err = max(abs(g_ex - g_n))
  # plt.plot(g_ex,'-')
  # plt.plot(g_n, 'r-')
  # plt.show(block=False)
  # plt.figure()
  return err

def normErrInf():
  xh = np.loadtxt('out/mesh_one_ellipse_reduced1_0_300_xh_tot.txt')
  yh = np.loadtxt('out/mesh_one_ellipse_reduced1_0_300_yh_tot.txt')
  return sg.Pointset(xh + 1j * yh)
