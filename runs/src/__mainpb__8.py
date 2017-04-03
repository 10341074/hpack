from inverseproblem import *

#########

h = 200.0
nsrc = 20

R = ()
a = 1e-16
reg = 1
regmet = 'tikh'
solver = 'lu'

theta = np.pi / 2
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
  sd = sg.Segment(nsd, f_inargs = (sh.ellipse, (0, rd, 2*rd)), quad='ps')
  #sd = sg.Segment(nsd, Z_args = (sh.dZ, sh.dZp, sh.dZpp, ()), quad='gp', aff=(2,1j))

  bd1 = sg.Boundary([sd])
  ld = sg.Layer([bd1], ly.layerpotSD)
  return (ld, so, sb)

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


plot.plot(x, y, ninv,'im')
ld.plot(p=True)
plt.show(block=False)

# R = np.array(R, float)
# Rreg = np.array(Rreg, float)

# plt.savefig('fig.pdf')
# plt.savefig('fig.png')
# plt.savefig('fig.ps')
# plt.savefig('fig.eps')
plt.savefig('fig_ninv.svg')
end = input('Press enter')

