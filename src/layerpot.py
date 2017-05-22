from scipy.linalg import circulant
import numpy as np

import __types__
import quadr

pi = np.pi
log = np.log
cos = np.cos
sin = np.sin
real = np.real
symmflagval = -999.; # all diag vals of this signifies symmetric - a hack
def scalar(a, b):
  return a.real * b.real + a.imag * b.imag
  # return np.real(a * np.conj(b))

def fundsol(r = [], k = 0):
  return -1.0 / 2 / pi * log(abs(r))

def fundsol_deriv(r = [], cosphi = 0 , k = 0): # without a -1
  return 1.0 / 2 / pi / r * cosphi

def phi(z0=0, z=[]):
  return - 1. / 2 / pi * log(abs(z - z0))
def phi_p(z0=0, z=[]):
  d = z - z0
  return - 1. / 2 / pi * d/abs(d)**2
def phi_n(z0=0, z=[], n=[]):
  return np.real(np.conj(n) * phi_p(z0, z))
def phi_theta(z0=0, z=[], theta=0):
  n = (np.cos(theta) + 1j * np.sin(theta) ) * np.ones((len(z)))
  return scalar(n, phi_p(z0, z))

def phi_x(z0=0, z=[]):
  return np.real(phi_p(z0, z))
def phi_y(z0=0, z=[]):
  return np.imag(phi_p(z0, z))

def phi_xx(z0=0, z=[]):
  return 1. / 2 / pi * 2 / abs(z-z0)**3 * np.real(z-z0) * np.real(z-z0) / abs(z-z0) - 1. / 2 / pi / abs(z-z0)**2
def phi_yy(z0=0, z=[]):
  return 1. / 2 / pi * 2 / abs(z-z0)**3 * np.imag(z-z0) * np.imag(z-z0) / abs(z-z0) - 1. / 2 / pi / abs(z-z0)**2
def phi_xy(z0=0, z=[]):
  return 1. / 2 / pi * 2 / abs(z-z0)**3 * np.real(z-z0) * np.imag(z-z0) / abs(z-z0)

def phi_hess(z0=0, z=[]):
  return np.array(
    [[ phi_xx(z0, z), phi_xy(z0, z)],
     [ phi_xy(z0, z), phi_yy(z0, z)]])

def phi_xxx(z0=0, z=[]):
  return 1. / 2 / pi * (-8) / abs(z-z0)**5 * np.real(z-z0) / abs(z-z0) * np.real(z-z0)**2 + 1. / 2 / pi * 4 / abs(z-z0)**4 * np.real(z-z0) + 1. / 2 / pi * 2/ abs(z-z0)**3 * np.real(z-z0) / abs(z-z0)
def phi_xxy(z0=0, z=[]):
  return 1. / 2 / pi * (-8) / abs(z-z0)**5 * np.imag(z-z0) / abs(z-z0) * np.real(z-z0)**2 + 1. / 2 / pi * 2/ abs(z-z0)**3 * np.imag(z-z0) / abs(z-z0)
def phi_yyx(z0=0, z=[]):
  return 1. / 2 / pi * (-8) / abs(z-z0)**5 * np.imag(z-z0) / abs(z-z0) * np.imag(z-z0)**2 + 1. / 2 / pi * 2/ abs(z-z0)**3 * np.real(z-z0) / abs(z-z0)
def phi_yyy(z0=0, z=[]):
  return 1. / 2 / pi * (-8) / abs(z-z0)**5 * np.imag(z-z0) / abs(z-z0) * np.imag(z-z0)**2 + 1. / 2 / pi * 4 / abs(z-z0)**4 * np.imag(z-z0) + 1. / 2 / pi * 2/ abs(z-z0)**3 * np.imag(z-z0) / abs(z-z0)

def phi_x_p(z0=0, z=[]):
  return phi_xx(z0, z) + 1j * phi_xy(z0, z)
def phi_y_p(z0=0, z=[]):
  return phi_xy(z0, z) + 1j * phi_yy(z0, z)

def phi_l(z0=0, z=[]):
  return phi_xx(z0, z) + phi_yy(z0, z)
def phi_l_p(z0=0, z=[]):
  return phi_xxx(z0, z) + phi_yyx(z0, z) + 1j * phi_xxy(z0, z) + 1j * phi_yyy(z0, z)

def phi_x_n(z0=0, z=[], n=[]):
  return np.real(np.conj(n) * phi_x_p(z0, z))
def phi_y_n(z0=0, z=[], n=[]):
  return np.real(np.conj(n) * phi_y_p(z0, z))
def phi_l_n(z0=0, z=[], n=[]):
  return scalar(phi_l_p(z0, z), n)

def circulant_T(a=[]):
  A = circulant(a)
  A = A.T
  return A

def layerpotS(k=0, s=[], t=[], o=[], nodiag=0):
  slf = 0
  if t == []:
    slf = 1
    t = s
  M = len(t.x)
  N = len(s.x)

  d = np.array([t.x for k in range(N)])
  d = d.T
  d = d - np.array([s.x for k in range(M)])
  r = abs(d)
  if slf or nodiag:
    r[np.diag_indices(N)] = symmflagval

  A = fundsol(r, k)
  sp = s.speed / 2 / pi # note: 2pi converts to speed wrt s in [0,2pi]
  if slf:
    S1 = -1. / 4 / pi # const M_1/2 of Kress w/out speed fac
    # A = A - S1 * circulant_T(log(4. * sin(pi / N * np.arange(N))**2 )) # A=D2=M_2/2
    A = A - S1 * circulant_T(np.concatenate(( [0], log(4. * sin(pi / N * np.arange(1,N))**2 ) )) ) # A=D2=M_2/2
    A[np.diag_indices(N)] = -log(sp) / 2 / pi # diag vals propto curvature

    A = S1 * circulant_T(quadr.kress_Rjn(float(N)/2)) + 2. * pi / N * A
    A = A.dot(np.diag(sp))
  else:
    A = A.dot(np.diag(s.w))
  return A

def layerpotD(k=0, s=[], t=[], o=[]):
  slf = 0
  if t == []:
    slf = 1
    t = s
  M = len(t.x)
  N = len(s.x)
  
  d = np.array([t.x for k in range(N)])
  d = d.T
  d = d - np.array([s.x for k in range(M)])
  r = abs(d)
  if slf:
    r[np.diag_indices(N)] = symmflagval

  n = np.array([s.nx for k in range(M)])
  cosphi = np.real(np.conj(n) * d) / r;
  
  A = fundsol_deriv(r, cosphi, k)

  sp = s.speed / 2 / pi
  if slf:
    A[np.diag_indices(N)] = -s.kappa / 4 / pi
    A = A.dot(np.diag(s.w))
  else:
    A = A.dot(np.diag(s.w))
  return A

def layerpotSD(k=0, s=[], t=[], o=[]):
  slf = 0
  if t == []:
    slf = 1
    t = s
  M = len(t.x);
  N = len(s.x);
  
  d = np.array([t.x for k in range(N)])
  d = d.T
  d = d - np.array([s.x for k in range(M)])
  r = abs(d)
  if slf:
    r[np.diag_indices(N)] = symmflagval

  n = np.array([t.nx for k in range(N)])
  n = n.T
  cosphi = - np.real(np.conj(n) * d) / r;
  
  A = fundsol_deriv(r, cosphi, k)

  sp = s.speed / 2 / pi
  if slf:
    A[np.diag_indices(N)] = -s.kappa / 4 / pi
    A = A.dot(np.diag(s.w))
  else:
    A = A.dot(np.diag(s.w))
  return A

def layerpotDD(k=0, s=[], t=[], o=[]):
  slf = 0
  if t == []:
    slf = 1
    t = s
    print('Warning: layerpotDD self not implemented')
  M = len(t.x)
  N = len(s.x)
  
  d = np.array([t.x for k in range(N)])
  d = d.T
  d = d - np.array([s.x for k in range(M)])
  r = abs(d)
  if slf:
    r[np.diag_indices(N)] = symmflagval

  ny = np.array([s.nx for k in range(M)])
  # cosphi = np.real(np.conj(n) * d) / r;
  
  nx = np.array([t.nx for k in range(N)])
  nx = nx.T

  A = -2. / r**3 * 1. / r * scalar(d, nx) * scalar(d, ny) + 1. / r**2 * scalar(ny, nx)
  A = 1. / 2 / pi * A

  sp = s.speed / 2 / pi
  if slf:
    A[np.diag_indices(N)] = -s.kappa / 4 / pi
    A = A.dot(np.diag(s.w))
  else:
    A = A.dot(np.diag(s.w))
  return A

    
def layerpotSDnow(k=0, s=(), t=()):
  M = len(t.x);
  N = len(s.x);
  
  d = np.array([t.x for k in range(N)])
  d = d.T
  d = d - np.array([s.x for k in range(M)])
  r = abs(d)

  n = np.array([t.nx for k in range(N)])
  n = n.T
  cosphi = - np.real(np.conj(n) * d) / r;
  
  A = fundsol_deriv(r, cosphi, k)
  return A

def layerpotSnow(k=0, s=(), t=()):
  M = len(t.x)
  N = len(s.x)

  d = np.array([t.x for k in range(N)])
  d = d.T
  d = d - np.array([s.x for k in range(M)])
  r = abs(d)

  A = fundsol(r, k)
  return A
