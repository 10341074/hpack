from __load__ import *
from inverseproblem import *
#import __mainpb__ as _m # for test

tt = time.time()
tc = time.clock()
#########

h = 200.0
nsrc = 20

R = ()
a = 1e-14
reg = 1
regmet = 'tikh'
solver = 'lu'

theta = 0
#########

def x3domain():
  ns = nsrc
  
  rb  = 10
  nsb = rb * ns
  sb = sg.Segment(nsb, f_inargs = (sh.circle, (0, rb)), quad='ps')

  ro  = 6
  nso = ro * ns
  so = sg.Segment(nso, f_inargs = (sh.circle, (0, ro)), quad='ps')

  rd  = 1
  nsd = rd * ns
  sd = sg.Segment(nsd, f_inargs = (sh.circle, (0, rd)), quad='ps')
  sd = sg.Segment(nsd, f_inargs = (sh.ellipse, (0, rd, 1.5*rd)), quad='ps')
  # sd = sg.Segment(nsd, Z_args = (sh.dZ, sh.dZp, sh.dZpp, ()), quad='gp', aff=(0, 0.8 + 0.8j))

  bd1 = sg.Boundary([sd])
  ld = sg.Layer([bd1], ly.layerpotSD)
  return (ld, so, sb)

def method_gap():
  c = 1.0 * (h + 1) / (h - 1)
  print(c)
  ld, so, sb = x3domain()
  nsd, nso, nsb = ld.n, so.n, sb.n

  allpsi = computeallpsi(ld, sb, c)
  # dpb.plotdpb(ld, sb.x[0], (-2, 2, 100), psi = allpsi[:,0], t='im')

  R, U, U_nu = computeR(allpsi, ld, so, sb)
  _gap = gap_init(R, a, reg, regmet, solver)

  # bb = sg.Boundary([sb])
  # lb = sg.Layer(b=[bb])
  x, y, pp = meshgrid((-2, 2 , 40))

  (ninv, res, nsolgap) = computeallsolsgap(_gap, pp, R, U, U_nu, so, theta)

  ninv01 = ninv / max(ninv)
 
  plot.plot(x, y, ninv, 'im')
  ld.plot(p=True)
  plt.show(block=False)

  # R = np.array(R, float)
  # Rreg = np.array(Rreg, float)

  # plt.savefig('fig.pdf')
  # plt.savefig('fig.png')
  # plt.savefig('fig.ps')
  # plt.savefig('fig.eps')
  plt.savefig('fig_ninv.svg')
  return ninv

def method_NtoD():
  c = 1.0 * (h + 1) / (h - 1)
  print(c)
  ld, so, sb = x3domain()
  nsd, nso, nsb = ld.n, so.n, sb.n
  print('theta = ', theta)

  LL0 = computeLL0(ld, so, sb, c, testset=0)
  _NtoD = NtoD_init(LL0, a, reg, regmet, solver)
  
  x, y, pp = meshgrid((-2, 2 , 40))
  
  (ninv, res, nsolgap) = computeallsolsNtoD(_NtoD, pp, LL0, so, theta)
  ninv01 = ninv / max(ninv)

  plot.plot(x, y, ninv, 'cf')
  ld.plot(p=True)
  plt.show(block=False)
  plt.savefig('fig_ninv.svg')
  return ninv

c = 1.0 * (h + 1) / (h - 1)
print(c)
ld, so, sb = x3domain()
nsd, nso, nsb = ld.n, so.n, sb.n

# method_gap()
thetav = np.pi * np.array([0], float)
for theta in np.array(thetav):
  ninv = method_NtoD()
tt = time.time() - tt
tc = time.clock() - tc

print('time wall-clock = ', tt)
print('time clock = ', tc)

end = input('Press enter')

def test_NtoD():
  g = ly.phi_n(-8, so.x, so.n)
  #g = np.ones(so.n)
  psi0 = dpb.mapNtoD0(so, g)
  bo = sg.Boundary([so])
  lo = sg.Layer([bo], dns=psi0)
  
  dpb.plotdpb(lo, (), (-8, 8, 100), t='im')

  psi = dpb.mapNtoD(so, ld, g, c)
  lo.dns = psi[0:lo.n]
  ld.dns = psi[lo.n:]

  dpb.plotdpb(lo, (), (-8, 8, 100), t='im', l2=ld)
  return
