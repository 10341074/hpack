import numpy as np

import shapes as sh
import segment as sg


def one_ellipse(nsd):
  return sg.Layer([sg.Boundary([sg.Segment(nsd, f_inargs = (sh.ellipse, (0, 2, 1)), quad='ps') ])])
def two_ellipse(nsd):
  return sg.Layer([
      sg.Boundary([sg.Segment(nsd, f_inargs = (sh.ellipse, (0, 2, 1)), quad='ps', aff=(-1-1j, 0.4 + 0.4j)) ]), # nsd/2
      sg.Boundary([sg.Segment(nsd, f_inargs = (sh.ellipse, (0, 2, 1)), quad='ps', aff=(0.5 +1j, 0.4 - 0.4j)) ]) # nsd/2
      ])
def three_ellipse(nsd):
  return sg.Layer([
      sg.Boundary([sg.Segment(nsd, f_inargs = (sh.ellipse, (0, 2, 1)), quad='ps', aff=(-1+0.2j, 0.3*np.exp(1j * np.pi/6))) ]),
      sg.Boundary([sg.Segment(nsd, f_inargs = (sh.ellipse, (0, 2, 1)), quad='ps', aff=(0.5 +1j, 0.2*np.exp(1j * np.pi/4))) ]),
      sg.Boundary([sg.Segment(nsd, f_inargs = (sh.ellipse, (0, 2, 1)), quad='ps', aff=(0.5 -1j, 0.4*np.exp(1j * np.pi/2))) ])
      ])
def one_kite(nsd):
  return sg.Layer([sg.Boundary([sg.Segment(nsd, Z_args = (sh.kZ, sh.kZp, sh.kZpp, ()), quad='ps', aff=(0, 0.8 + 0.8j)) ])])
def one_drop(nsd):
  return sg.Layer([sg.Boundary([sg.Segment(nsd, Z_args = (sh.dZ, sh.dZp, sh.dZpp, ()), quad='gp', aff=(0, 0.8 + 0.8j)) ])])
def one_triangle(nsd):
  return sg.Layer([sg.Boundary([ sg.poly([2, 1 +1.5j, -1 +1j], int(nsd/3)) ])])
#######################
def sg_one_ellipse(nsd):
  return sg.Segment(nsd, f_inargs = (sh.ellipse, (0, 2, 1)), quad='ps')
def sg_one_kite(nsd):
  return sg.Segment(nsd, Z_args = (sh.kZ, sh.kZp, sh.kZpp, ()), quad='ps', aff=(0, 0.8 + 0.8j))
def sg_one_drop(nsd):
  return sg.Segment(nsd, Z_args = (sh.dZ, sh.dZp, sh.dZpp, ()), quad='gp', aff=(0, 0.8 + 0.8j))
def sg_one_triangle(nsd):
  return sg.poly([2, 1 +1.5j, -1 +1j], nsd) # nsd/3
# def two_ellipse(nsd):
#   return sg.Layer([
#       sg.Boundary([sg.Segment(int(nsd/2), f_inargs = (sh.ellipse, (0, 2, 1)), quad='ps', aff=(-1-1j, 0.4 + 0.4j)) ]),
#       sg.Boundary([sg.Segment(int(nsd/2), f_inargs = (sh.ellipse, (0, 2, 1)), quad='ps', aff=(0.5 +1j, 0.4 - 0.4j)) ])
#       ])
# def one_kite(nsd):
#   return sg.Layer([sg.Boundary([sg.Segment(nsd, Z_args = (sh.kZ, sh.kZp, sh.kZpp, ()), quad='ps', aff=(0, 0.8 + 0.8j)) ])])

# def one_drop(nsd):
#   return sg.Layer([sg.Boundary([sg.Segment(nsd, Z_args = (sh.dZ, sh.dZp, sh.dZpp, ()), quad='gp', aff=(0, 0.8 + 0.8j)) ])])
# def one_triangle(nsd):
#   return sg.Layer([sg.Boundary([ sg.poly([2, 1 +1.5j, -1 +1j], int(nsd/3)) ])])



def example(index = 'one_ellipse', nsb = 80, nso = 100, nsd = 40, rb = 5, ro = 3):
  sb = sg.Segment(nsb, f_inargs = (sh.circle, (0, rb)), quad='ps')
  so = sg.Segment(nso, f_inargs = (sh.circle, (0, ro)), quad='ps')
  switcher_d = {
    'one_ellipse' : one_ellipse,
    'two_ellipse' : two_ellipse,
    'one_kite' : one_kite,
    'one_drop' : one_drop,
    'one_triangle' : one_triangle
  }
  return (sb, so, switcher_d[index](nsd))

def domain(nsb=80, nso=80, nsd=40, rb=10, ro=3):
  ns = nsrc
  
  rb  = 10
  nsb = rb * ns
  if nsb == 0:
    nsb = 160
  sb = sg.Segment(nsb, f_inargs = (sh.circle, (0, rb)), quad='ps')

  ro  = 3
  nso = ro * ns
  if nso == 0:
    nso = 180
  so = sg.Segment(nso, f_inargs = (sh.circle, (0, ro)), quad='ps')
  # so = sg.Segment(nso, f_inargs = (sh.ellipse, (0, 4, 3)), quad='ps')

  rd  = 1
  if nsd == 0:
    nsd = 80
  # sd = sg.Segment(nsd, f_inargs = (sh.circle, (0, rd)), quad='ps')
  sd = sg.Segment(nsd, f_inargs = (sh.ellipse, (0, 2*rd, rd)), quad='ps')
  # sd = sg.Segment(nsd, Z_args = (sh.dZ, sh.dZp, sh.dZpp, ()), quad='gp', aff=(0, 0.8 + 0.8j))
  # sd = sg.Segment(nsd, Z_args = (sh.kZ, sh.kZp, sh.kZpp, ()), quad='ps', aff=(0, 0.8 + 0.8j))
  sd1 = sg.Segment(int(nsd/2), f_inargs = (sh.ellipse, (0, 2*rd, rd)), quad='ps', aff=(-1-1j, 0.4 + 0.4j))
  sd2 = sg.Segment(int(nsd/2), f_inargs = (sh.ellipse, (0, 2*rd, rd)), quad='ps', aff=(0.5 +1j, 0.4 - 0.4j))
  
  # sd1 = sg.Segment(nsd, f_inargs = (sh.circle, (0, rd)), quad='ps', sign=-1)
  # sd2 = sg.Segment(nsd, f_inargs = (sh.circle, (0, 2*rd)), quad='ps')
  
  bd1 = sg.Boundary([sd1])
  bd2 = sg.Boundary([sd2])

  b = sg.Boundary([sd])
  # bd2 = sg.poly([2, 1 +1.5j, -1 +1j], n=20)
  ld = sg.Layer([bd1, bd2], ly.layerpotSD)
  # ld = sg.Layer([b], ly.layerpotSD)
  return (ld, so, sb)

