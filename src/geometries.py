import numpy as np

import shapes as sh
import segment as sg
# -----------------------------------------------------------------------------------------------
def lay_one_ellipse(nsd):
  return sg.Layer([sg.Boundary([sg.Segment(nsd, f_inargs = (sh.ellipse, (0, 2, 1)), quad='ps') ])])
def lay_two_ellipse(nsd):
  return sg.Layer([
      sg.Boundary([sg.Segment(nsd, f_inargs = (sh.ellipse, (0, 2, 1)), quad='ps', aff=(-1-1j, 0.4 + 0.4j)) ]), # nsd/2
      sg.Boundary([sg.Segment(nsd, f_inargs = (sh.ellipse, (0, 2, 1)), quad='ps', aff=(0.5 +1j, 0.4 - 0.4j)) ]) # nsd/2
      ])
def lay_three_ellipse(nsd):
  return sg.Layer([
      sg.Boundary([sg.Segment(nsd, f_inargs = (sh.ellipse, (0, 2, 1)), quad='ps', aff=(-1+0.2j, 0.3*np.exp(1j * np.pi/6))) ]),
      sg.Boundary([sg.Segment(nsd, f_inargs = (sh.ellipse, (0, 2, 1)), quad='ps', aff=(0.5 +1j, 0.2*np.exp(1j * np.pi/4))) ]),
      sg.Boundary([sg.Segment(nsd, f_inargs = (sh.ellipse, (0, 2, 1)), quad='ps', aff=(0.5 -1j, 0.4*np.exp(1j * np.pi/2))) ])
      ])
def lay_one_kite(nsd):
  return sg.Layer([sg.Boundary([sg.Segment(nsd, Z_args = (sh.kZ, sh.kZp, sh.kZpp, ()), quad='ps', aff=(0, 0.8 + 0.8j)) ])])
def lay_one_drop(nsd):
  return sg.Layer([sg.Boundary([sg.Segment(nsd, Z_args = (sh.dZ, sh.dZp, sh.dZpp, ()), quad='gp', aff=(0, 0.8 + 0.8j)) ])])
def lay_one_triangle(nsd):
  return sg.Layer([sg.Boundary([ sg.poly([2, 1 +1.5j, -1 +1j], int(nsd/3)) ])])
# -----------------------------------------------------------------------------------------------
def sg_one_circle(ns):
  return sg.Segment(ns, f_inargs = (sh.circle, (0, 2)), quad='ps')
def sg_one_ellipse(ns):
  return sg.Segment(ns, f_inargs = (sh.ellipse, (0, 2, 1)), quad='ps')
def sg_one_kite(ns):
  return sg.Segment(ns, Z_args = (sh.kZ, sh.kZp, sh.kZpp, ()), quad='ps', aff=(0, 0.8 + 0.8j))
def sg_one_drop(ns):
  return sg.Segment(ns, Z_args = (sh.dZ, sh.dZp, sh.dZpp, ()), quad='gp', aff=(0, 0.8 + 0.8j))
def sg_one_triangle(ns):
  return sg.poly([2, 1 +1.5j, -1 +1j], ns) # ns/3
# -----------------------------------------------------------------------------------------------
def bod_template_c(name_sd = 'one_ellipse', nsb = 80, nso = 100, nsd = 40, rb = 5, ro = 3, ns=0, quad='ps'):
  '''
  This function returns a tuple (sb, so, sd) with three segment.Layer objects of size 'nsb', 'nso', 'nsd'
  sb, so are circles, sd can be set among some examples with 'name_sd'
  ('ns' can be used to scale nsb, nso proportionally to radii 'rb', 'ro')
  '''
  nsb, nso = (nsb, nso) if ns <= 0 else (ns * rb, ns * ro)
  sb = sg.Segment(nsb, f_inargs = (sh.circle, (0, rb)), quad=quad)
  so = sg.Segment(nso, f_inargs = (sh.circle, (0, ro)), quad=quad)
  switcher_sd = {
    'one_ellipse' : lay_one_ellipse,
    'two_ellipse' : lay_two_ellipse,
    'three_ellipse' : lay_three_ellipse,
    'one_kite' : lay_one_kite,
    'one_drop' : lay_one_drop,
    'one_triangle' : lay_one_triangle
  }
  return (sb, so, switcher_sd[name_sd](nsd))

example = bod_template_c

