import scipy.linalg as linalg
import numpy as np

def _solve(_m, _rhs):
  return linalg.solve(_m, _rhs)
def _lu(_m, _rhs):
  return linalg.lu_solve(_m, _rhs)
def _lstsq(_m, _rhs):
  return linalg.lstsq(_m, _rhs)[0]

class Solver:
  def __init__(self, A=(), RHS_f=(), RHS_fcom=(), RHS_args=(), BX=(), BY=()):
    self.A = A
    self.RHS_f = RHS_f
    self.RHS_fcom = RHS_fcom
    self.RHS_args = RHS_args

    self.solver_a = ()
    self.solver_f = ()

    if BX == ():
      self.BX = np.eye(A.shape[1])
    else:
      self.BX = BX
    if BY == ():
      self.BY = np.eye(A.shape[0])
    else:
      self.BY = BY
