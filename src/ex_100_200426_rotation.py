from src import *
'''
  I want to analyze the effect of a small rotation introduced by the shift in the quadrature rule

'''
# ----------------
radiusb = 4
nb, no, nd = 300, 300, 200
# ----------------

def init(quad):
  p = m.EIT()
  p.domain('one_ellipse', nso = no, quad=quad)
  p.meshgrid((-3, 3, 80))
  # -----------------------------------------
  p.alpha = 1e-15
  p.noiselevel = 0
  # p.c = 100
  return p

p1 = init('p')
p2 = init('ps')

p1.ld = sg.Layer([sg.Boundary([sg.Segment(nd, Z_args = (sh.kZ, sh.kZp, sh.kZpp, ()), quad='p') ])])
p2.ld = sg.Layer([sg.Boundary([sg.Segment(nd, Z_args = (sh.kZ, sh.kZp, sh.kZpp, ()), quad='p', aff=(0, np.exp(1j * np.pi / 300))) ])])
# p1.ld = sg.Layer([sg.Boundary([sg.Segment(nd, f_inargs = (sh.ellipse, (0, 2, 1)), quad='p') ])])
# p2.ld = sg.Layer([sg.Boundary([sg.Segment(nd, f_inargs = (sh.ellipse, (0, 2, 1)), quad='p', aff=(0, np.exp(1j * np.pi / 300))) ])])

p1.so = sg.Segment(nd, f_inargs = (sh.ellipse, (0, 3, 3.2)), quad='p')
p2.so = sg.Segment(nd, f_inargs = (sh.ellipse, (0, 3, 3.2)), quad='p', aff=(0, np.exp(1j * np.pi / 300)))

p1.sb = sg.Segment(nb, f_inargs = (sh.circle, (0, radiusb)), quad='p')
p2.sb = sg.Segment(nb, f_inargs = (sh.circle, (0, radiusb)), quad='ps')

p1.theta = 0
p2.theta = np.pi / 300


def rg(p):
  p.rg_solver()
  p.ipb()
  fig = plt.figure()
  plt.contour(p.x, p.y, p.z.reshape(p.y.size, p.x.size), 30)
  p.ld.plot(lw = 1.2, ms = 1)
  p.so.plot(lw = 1.2, ms = 1)
  plt.axis('equal')
  plt.axis('square')
  plt.show(block=False)

rg(p1)
rg(p2)

input('press')