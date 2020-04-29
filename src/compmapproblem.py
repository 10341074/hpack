import matplotlib.pyplot as plt
import numpy as np
from numpy.linalg import norm
from numpy.linalg import cond
import scipy.linalg as linalg
import numpy.linalg
import scipy.stats
import numpy.random
import time

from __types__ import *
import shapes as sh
import layerpot as ly
import segment as sg
import plot
import directproblem as dpb

import linfunc as linf
import lintype as lint
import setups

import debug_globals as dbgg
from debug import printname
verbose = 1

# -------------------------------------------------------------
def printname_compL(func):
  return printname(func, "compL")
@printname_compL
def computeL0(so, T):
  if T == ():
    T = so.BX
  allpsi = dpb.mapNtoD0(so, T, so.s0)
  # ---------------------------------------  
  Lo = ly.layerpotS(s=so)
  L0 = Lo.dot(allpsi)
  means = np.ones(so.n).dot(np.diagflat(so.w).dot(L0)) / sum(so.w)
  L0 = L0 - np.array([means for k in range(so.n)])
  return L0
@printname_compL
def computeL0B(so, T, L0=()):
  if L0 == ():
    print('  not passed L0')
    L0 = computeL0(so, T)
  else:
    print('  passed L0')
  # return so.BY.T.dot(np.diagflat(so.w).dot(L0))
  return so.BYinv.dot(L0)
@printname_compL
def computeL(ld, so, T, c):
  if T == ():
    T = so.BX
  allpsi = dpb.mapNtoD(so, ld, T, c, so.s0)
  # ---------------------------------------  
  Lo = ly.layerpotS(s=so)
  Ld = ly.layerpotS(s=ld, t=so)
  L = Lo.dot(allpsi[0:so.n]) + Ld.dot(allpsi[so.n::])
  # means = sum(np.diagflat(so.w).dot(L)) / sum(so.w) # correct? strange sum by rows
  means = np.ones(so.n).dot(np.diagflat(so.w).dot(L)) / sum(so.w)
  L = L - np.array([means for k in range(so.n)])
  return L
@printname_compL
def computeLB(ld, so, T, c, L=()):
  if L == ():
    L = computeL(ld, so, T, c)
    print('  not passed L')
  else:
    print('  passed L')
  # return so.BY.T.dot(np.diagflat(so.w).dot(L))
  return so.BYinv.dot(L)
@printname_compL
def computeLL0(ld, so, T, c, L0=(), L=()):
  if T == ():
    T = so.BX
  # Kp = ly.layerpotSD(s=so)
  # (nu, s0) = linf.eigmaxpowerw(A=Kp, s=so)
  # if testset == 0:
  #   S = linf.gramschmidtw(so, s0=s0)
  #   T = S[:, 1::]
  # if testset == 1 or testset == 3:
  #   V1_nu = ly.layerpotSD(s=sb, t=so)
  #   T = V1_nu
  # if testset == 2 or testset == 3:
  #   V2_nu = ly.layerpotDD(s=sb, t=so)
  #   T = V2_nu
  # if testset == 3:
  #   T = np.concatenate((V1_nu.T, V2_nu.T)).T
  if L0 == ():
    print('  not passed L0')
    L0 = computeL0(so, T)
  else:
    print('  passed L0')
  if L == ():
    print('  not passed L')
    L = computeL(ld, so, T, c)
  else:
    print('  passed L')
  LL0 = L - L0
  return LL0
@printname_compL
def computeLL0B(ld, so, T, c, L0B=(), LB=()):
  if L0B == ():
    L0B = computeL0B(so, T)
    print('computing L0B/LL0B')
  else:
    print('  taking L0B/LL0B from L0B')
  if LB == ():
    LB = computeLB(ld, so, T, c)
    print('computing LB/LL0B')
  else:
    print('  taking LB/LL0B from LB')

  LL0 = LB - L0B
  return LL0
@printname_compL
def computeLLdiff(ld, so, T, c):
  if T == ():
    T = so.BX
  allpsi0 = dpb.mapNtoD0(so, T, so.s0) 
  Lo = ly.layerpotSD(s=so, t =ld)
  rhsdiff = Lo.dot(allpsi0)
  allpsi = dpb.mapNtoDdiff(so, ld, rhsdiff, c, so.s0)
  Lo = ly.layerpotS(s=so)
  Ld = ly.layerpotS(s=ld, t=so)
  L = Lo.dot(allpsi[0:so.n]) + Ld.dot(allpsi[so.n::])
  # means = sum(np.diagflat(so.w).dot(L)) / sum(so.w) # correct? strange sum by rows
  means = np.ones(so.n).dot(np.diagflat(so.w).dot(L)) / sum(so.w)
  L = L - np.array([means for k in range(so.n)])
  return L
@printname_compL
def computeLLBdiff(ld, so, T, c):
  if T == ():
    T = so.B[:, 1:]
  L = computeLLdiff(ld, so, T, c)
  return so.BY.T.dot(np.diagflat(so.w).dot(L))
# ----------------------------------------------------------------------------------
