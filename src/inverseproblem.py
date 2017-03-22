# from domain import *
import matplotlib.pyplot as plt
import numpy as np
from scipy.linalg import solve
from numpy.linalg import norm
from numpy.linalg import cond


import my_type
import layerpot as ly
import segment as sg
import plot
import directproblem as dpb

def computeRreg(R, a, nsb):
  return a * np.eye(nsb) + R.T.dot(R)

def computeRHS(U, U_nu, z0):
  F = ly.phi_x(z0, so.x)
  F_nu = ly.phi_x_n(z0, so.x, so.nx)
  return U.dot(F_nu) - U_nu.dot(F)

def computeRHSreg(R, U, U_nu, z0, RHS=[]):
  if RHS == []:
    RHS = computeRHS(U, U_nu, z0)
  return R.T.dot(RHS)

h = 20
n = 20
c = 1. * (h + 1)/(h-1)

ns = n

rb  = 10
nsb = rb * ns
sb = sg.Segment(nsb, f=sg.circle, inargs=(1, rb), periodic=True)

ro  = 6
nso = ro * ns
so = sg.Segment(nso, f=sg.circle, inargs=(1j, ro), periodic=True)

rd  = 1
nsd = rd * ns
sd = sg.Segment(nsd, f=sg.circle, inargs=(0, rd), periodic=True)
#sd = sg.Segment(nsd, f=sg.ellipse, inargs=(0, rd, 1.5*rd), periodic=True)

# (K' + 0.5 * c(h) * I) psi = -phi_nu

# s = sg.Segment(n, f = sg.circle, inargs = (0, 1))
# b = sg.Boundary([s])
# p = sg.Pointset(s.x)

ld = sg.Layer(sd, ly.layerpotSD)
A = dpb.directA(ld, c)

allpsi = np.empty((nsd, nsb), float)
for k in range(nsb):
  z0 = sb.x[k]
  rhs = dpb.directrhs(ld.b, z0)
  allpsi[:, k] = dpb.directpb(z0, ld, c, A, rhs)

# Kp = l.exp(0, l.b, [], [])
# A = Kp
# A[np.diag_indices(len(s.x))]=A[np.diag_indices(len(s.x))] + 0.5 * c

# z0 = 4
# d = s.x - z0
# r = abs(d)
# cosphi = np.real(np.conj(s.nx) * d) / r
# rhs = - ly.fundsol_deriv(r, cosphi, 0)

# psi = solve(A, rhs)


z0 = sb.x[0]
x, y, pp = plot.meshgrid((-2, 2, 20))
G = sg.Layer(b=sd, exp=ly.layerpotS, dns=allpsi[:,0])

pp = sg.Pointset(pp)

vv1 = sg.eval_layer(G, pp)
vv2 = ly.fundsol(abs(pp.x - z0), 0)
h = plt.figure()
plot.plot(h, x, y, vv2 + vv1)
sd.plot(p=True)
# so.plot(p=True)
# sb.plot(p=True)
plt.show(block=True)


U_psi = np.empty((nso, nsb), float)
U_psi_nu = np.empty((nso, nsb), float)
# kerS = sg.eval_layer(ld, so, exp = ly.layerpotS)
# kerSD = sg.eval_layer(ld, so, exp = ly.layerpotSD)
kerS = ly.layerpotS(0, sd, so)
kerSD = ly.layerpotSD(0, sd, so)

for k in range(nsb):
  U_psi[:,k] = kerS.dot(allpsi[:,k])
  U_psi_nu[:,k] = kerSD.dot(allpsi[:,k])
# U_psi = np.diag(so.w) * U_psi
# U_psi_nu = np.diag(so.w) * U_psi_nu
U_psi = U_psi.T.dot(np.diag(so.w))
U_psi_nu = U_psi_nu.T.dot(np.diag(so.w))

lo = sg.Layer(so, ly.layerpotS)
U = ly.layerpotS(0, so, sb)
U_nu = ly.layerpotD(0, so, sb)

# U = U + U_psi.';
# U_nu = U_nu + U_psi_nu.';

U = U + U_psi
U_nu = U_nu + U_psi_nu

lb = sg.Layer(sd, ly.layerpotS)
V = ly.layerpotS(0, sb, so)
V_nu = ly.layerpotSD(0, sb, so)

R = U.dot(V_nu) - U_nu.dot(V)


normf = np.empty(len(pp.x), float)
for k in range(len(pp.x)):
  z0 = pp.x[k]
  #F_stored(:,k) = F;
  #F_nu_stored(:,k) = F_nu;
  #RHS_stored(:,k) = RHS;

  #zeta = A \ RHS;
  a = 1e-14;

  Rreg = computeRreg(R, a, nsb)
  RHS = computeRHS(U, U_nu, z0)
  RHSreg = computeRHSreg(R, U, U_nu, z0, RHS)
    
  zeta = solve(Rreg, RHSreg)
  #sum(F_nu .* so.w.'); % check null
    
  normf[k] = norm(RHS) / norm(zeta)
  #normf(k) = real(z0) * imag(z0) + real(z0);

h = plt.figure()
plot.plot(h, x, y, normf)
sd.plot(p=True)
plt.show(block=True)

R = np.array(R, float)
Rreg = np.array(Rreg, float)


