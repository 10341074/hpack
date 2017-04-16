import scipy.linalg as linalg

def _solve(_m, _rhs):
  return linalg.solve(_m, _rhs)
def _lu(_m, _rhs):
  return linalg.lu_solve(_m, _rhs)
def _lstsq(_m, _rhs):
  return linalg.lstsq(_m, _rhs)[0]

class Solver:
  def __init__(self, A=(), RHS_f=(), RHS_args=()):
    self.A = A
    self.RHS_f = RHS_f
    self.RHS_args = RHS_args

    self.solver_a = ()
    self.solver_f = ()
