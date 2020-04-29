import scipy.linalg as linalg
import numpy as np

def _solve(_m, _rhs):
  return linalg.solve(_m, _rhs)
def _lu(_m, _rhs):
  return linalg.lu_solve(_m, _rhs)
def _lstsq(_m, _rhs):
  return linalg.lstsq(_m, _rhs)[0]


class SolverSavings:
  def __init__(self, num_testpoints, Aregularized_shape):
    self.zeta    = np.zeros((0, num_testpoints, Aregularized_shape[1]), float)
    self.sol     = np.zeros((0, num_testpoints, Aregularized_shape[1]), float)
    self.alpha   = np.zeros((0, num_testpoints                       ), float)
    self.disc    = np.zeros((0, num_testpoints                       ), float)
    self.disc_p  = np.zeros((0, num_testpoints                       ), float)
    self.ratio   = np.zeros((0, num_testpoints                       ), float)
    self.alpha_l = np.zeros((   num_testpoints), float)
    self.alpha_r = np.zeros((   num_testpoints), float)
    self.rhs     = np.zeros((0, num_testpoints, Aregularized_shape[0]), float)
  def add_iterations_alpha(self, it_alpha, num_testpoints, Aregularized_shape):
    self.zeta    = np.append(self.zeta,   np.zeros((it_alpha, num_testpoints, Aregularized_shape[1]), float), axis=0)
    self.sol     = np.append(self.sol,    np.zeros((it_alpha, num_testpoints, Aregularized_shape[1]), float), axis=0)
    self.alpha   = np.append(self.alpha,  np.zeros((it_alpha, num_testpoints                       ), float), axis=0)
    self.disc    = np.append(self.disc,   np.zeros((it_alpha, num_testpoints                       ), float), axis=0)
    self.disc_p  = np.append(self.disc_p, np.zeros((it_alpha, num_testpoints                       ), float), axis=0)
    self.ratio   = np.append(self.ratio,  np.zeros((it_alpha, num_testpoints                       ), float), axis=0)
    self.rhs     = np.append(self.rhs,    np.zeros((it_alpha, num_testpoints, Aregularized_shape[0]), float), axis=0)


class Solver:
  def __init__(self, A, Areg_func, RHS_f, RHS_fcom, RHS_args, BX, BY, alpha, delta, w, testpoints, solvertype='s'):
    self.A = A
    self.Ar = Areg_func(A=A, alpha=alpha)
    
    ''' current alpha: it's updated during computations '''
    self.alpha = alpha
    self.w = w
    
    self.RHS_f = RHS_f # function
    self.RHS_fcom = RHS_fcom
    self.RHS_args = RHS_args

    self.solver_a = ()
    self.solver_f = ()

    self.BX = np.eye(A.shape[1]) if BX == () else BX
    self.BY = np.eye(A.shape[0]) if BY == () else BY

    ''' delta is used for noise '''
    self.delta = delta
    ''' original alpha: saved to restart morozov iterations with different test point '''
    self.alpha_orig   = alpha
    self.alpha_l_orig = 1e-30
    self.alpha_r_orig = alpha
    ''' savings '''
    self.save = SolverSavings(len(testpoints.x), self.Ar.shape)
    ''' init arrays '''
    self.init_solver_arrays(solvertype)
    return
  def init_solver_arrays(self, solvertype):
    self.solvertype = solvertype
    if solvertype == 's':
      self.solver_a = self.Ar
      self.solver_f = _solve
    elif solvertype == 'lu':
      (lu, piv) = linalg.lu_factor(self.Ar)
      self.solver_a = (lu, piv)
      self.solver_f = _lu
    elif solvertype == 'lstsq':
      self.solver_a = self.Ar
      self.solver_f = _lstsq # to change
    else:
      raise Exception('Error: not valid solvertype')
    return
  def add_iterations_alpha(self, it_alpha, num_testpoints):
    self.save.add_iterations_alpha(it_alpha, num_testpoints, self.Ar.shape)
  