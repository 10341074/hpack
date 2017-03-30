# from domain import *
import matplotlib.pyplot as plt
import numpy as np
from scipy.linalg import solve
from scipy.linalg import lstsq
from numpy.linalg import norm
from numpy.linalg import cond
import scipy.linalg as linalg

from __types__ import *
import layerpot as ly
import segment as sg
import plot
import directproblem as dpb

verbose=1

def computeRreg(R, a, nsb):
  return a * np.eye(R.shape[1]) + R.T.dot(R)

def computeRHS(U, U_nu, z0):
  F = ly.phi_l(z0, so.x)
  F_nu = ly.phi_l_n(z0, so.x, so.nx)
  return U.dot(F_nu) - U_nu.dot(F)

def computeRHSreg(R, U, U_nu, z0, RHS=[]):
  if RHS == []:
    RHS = computeRHS(U, U_nu, z0)
  return R.T.dot(RHS)

def computek(zeta, z0, ld, lb):
  F = ly.phi(z0, ld.x)
  F_nu = ly.phi_n(z0, ld.x, ld.nx)
  APsi = dpb.mapNtoD(ld, F_nu)
  Psi = F
  vz = ly.layerpotS(s=lb, t=ld)
  vz = vz.dot(zeta)
  i_v = (vz**2).dot(ld.w)
  i_A = (APsi**2).dot(ld.w)
  i_P = (Psi**2).dot(ld.w)
  i_AP = (APsi * Psi).dot(ld.w)
  # (AP - v + (AP^2 - 2*v*AP + A*v + P*v - A*P)^(1/2))/(P - v)
  # k = (i_AP - i_v + (i_AP**2 - 2*i_v*i_AP + i_A*i_v + i_P*i_v - i_A*i_P)**(0.5))/(i_P - i_v)
  nn = (1 - h)*vz - APsi + h * Psi
  nn = (nn**2).dot(ld.w)
  return nn

def domain():
  n = 40
  ns = n

  rb  = 10
  nsb = rb * ns
  sb = sg.Segment(nsb, f=sg.circle, inargs=(0, rb), quad='ps')

  ro  = 6
  nso = ro * ns
  so = sg.Segment(nso, f=sg.circle, inargs=(1, ro), quad='ps')

  rd  = 1
  nsd = rd * ns
  sd = sg.Segment(nsd, f=sg.circle, inargs=(0, rd), periodic=True)
  # sd = sg.Segment(nsd, f=sg.ellipse, inargs=(0.3, rd, 1.5*rd), periodic=True)
  # nsd = 60
  # sd = sg.Segment(nsd, Z=sg.dZ, Zp=sg.dZp, Zpp=sg.dZpp, args=[], quad='gp')
  sd = sg.Segment(nsd, f=sg.circle, inargs=(0, rd), quad='ps')

  bd1 = sg.Boundary([sd])
  ld = sg.Layer([bd1], ly.layerpotSD)
  return (ld, so, sb)

h = 200.0
c = 1. * (h + 1)/(h-1)
ld, so, sb = domain()
nsd, nso, nsb = ld.n, so.n, sb.n

def computeallpsi(ld, sb, c):
  nsd = ld.n
  nsb = sb.n
  allpsi = np.empty((nsd, nsb), float)
  A = dpb.directKpc(ld, c)
  for k in range(nsb):
    z0 = sb.x[k]
    rhs = dpb.directrhs(l=ld, z0=z0)
    allpsi[:, k] = dpb.directpb(l=ld, c=c, z0=z0, A_f = (A, linalg.solve), rhs=rhs)
  return allpsi

  
allpsi = computeallpsi(ld, sb, c)
dpb.plotdpb(ld, sb.x[0], (-2, 2, 100), psi = allpsi[:,0], t='im')


sd = ld.b[0]
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

U = ly.layerpotS(0, so, sb)
U_nu = ly.layerpotD(0, so, sb)

# U = U + U_psi.';
# U_nu = U_nu + U_psi_nu.';

U = U + U_psi
U_nu = U_nu + U_psi_nu

#lb = sg.Layer([sd], ly.layerpotS)
testset = 3
if testset == 1 or testset == 3:
  V1 = ly.layerpotS(0, sb, so)
  V1_nu = ly.layerpotSD(0, sb, so)
  R1 = U.dot(V1_nu) - U_nu.dot(V1)
  R = R1
if testset == 2 or testset == 3:
  V2 = ly.layerpotD(0, sb, so)
  V2_nu = ly.layerpotDD(0, sb, so)
  R2 = U.dot(V2_nu) - U_nu.dot(V2)
  R = R2
if testset == 3:
  R = np.concatenate((R1.T, R2.T)).T
if verbose:
  print('R shape ', R.shape)


bb = sg.Boundary([sb])
lb = sg.Layer(b=[bb])
x, y, pp = plot.meshgrid((-2, 2 , 20))
pp = sg.Pointset(pp)

normf = np.empty(len(pp.x), float)
nn = np.empty(len(pp.x), float)
a = 1e-12

Rreg = computeRreg(R, a, nsb)
(lu, piv) = linalg.lu_factor(Rreg)
for k in range(len(pp.x)):
  z0 = pp.x[k]
  #F_stored(:,k) = F;
  #F_nu_stored(:,k) = F_nu;
  #RHS_stored(:,k) = RHS;

  #zeta = A \ RHS;
  RHS = computeRHS(U, U_nu, z0)
  RHSreg = computeRHSreg(R, U, U_nu, z0, RHS)
  
  #zeta = lstq(Rreg, RHSreg)[0]
  #zeta = solve(Rreg, RHSreg)
  zeta = linalg.lu_solve((lu, piv), RHSreg)
  #sum(F_nu .* so.w.'); % check null

  
  normf[k] = norm(RHS) / norm(zeta)
  # if abs(z0) < 1e-10:
  #   print('k= ',computek(zeta, z0, ld2, lb))
  #nn[k] = computek(zeta, z0, ld, lb)
  #normf(k) = real(z0) * imag(z0) + real(z0);
  #normf[k] = max(abs(solve(Rreg, RHSreg) - zeta))


plot.plot(x, y, normf,'cf')
sd.plot(p=True)
#so.plot(p=True)
plt.show(block=True)

R = np.array(R, float)
Rreg = np.array(Rreg, float)

