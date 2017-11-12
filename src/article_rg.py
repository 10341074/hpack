from src import *
import test_plots

savefig = False
plotstats = False
#######################################################################################################################
no, nd = 120, 100
###########################################################################################################################
ellipses = [
    sg.Layer([sg.Boundary([sg.Segment(nd, f_inargs = (sh.ellipse, (0, 2, 1)), quad='ps') ])]), 
    sg.Layer([sg.Boundary([sg.Segment(nd, f_inargs = (sh.ellipse, (0, 1.8, 1.2)), quad='ps', aff=(1, np.exp(1j * np.pi/4))) ])]),
    sg.Layer([sg.Boundary([sg.Segment(nd, f_inargs = (sh.ellipse, (0, 0.1, 1.2)), quad='ps', aff=(-1, np.exp(1j * np.pi * 2/3))) ])])
]

def first_compare(ld = (), prename = (), name = () ):
  if ld == ():
    n = len(ellipses)
    prename, name = 'one_ellipse', 'ellipse'
  else:
    n = 1
  for k in range(0, n):
    p = m.EIT()
    p.domain('one_ellipse', nso = no)
    # p.ld = sg.Layer([sg.Boundary([sg.Segment(40, f_inargs = (sh.ellipse, (0, 2, 1)), quad='ps') ])])
    # p.ld = sg.Layer([sg.Boundary([sg.Segment(80, f_inargs = (sh.ellipse, (0, 1.5, 1)), quad='ps', aff=(1, 1.2*np.exp(1j * np.pi/4))) ])])
    if ld == ():
      p.ld = ellipses[k]
    else:
      p.ld = ld
    ###################################################################################################################
    p.sb = sg.Segment(120, f_inargs = (sh.circle, (0, 4)), quad='ps')
    p.meshgrid((-3, 3, 80))
    ############################################
    p.noiselevel = 0
    p.alpha = 1e-10
    p.rg_solver()
    p.rg_ipb()
    #############################################
    fig = plt.figure()
    # plot
    plt.contour(p.x, p.y, p.z.reshape(p.y.size, p.x.size), 30)
    p.ld.plot(lw = 1.2, ms = 1)
    p.so.plot(lw = 1.2, ms = 1)
    plt.axis('equal')
    plt.axis('square')
    plt.show(block=False)
    if savefig:
      plt.savefig('runs/fig-thesis/%s_rg_%s%s_alpha%s_no%snd%snb%s_radius%s.eps' %(prename, name, k, p.alpha, no, nd, p.sb.n, 4), bbox_inches='tight')
    ###################################################################################################################
    p.sb = sg.Segment(120, f_inargs = (sh.circle, (0, 6)), quad='ps')
    p.meshgrid((-3, 3, 80))
    ############################################
    p.noiselevel = 0
    p.alpha = 1e-10
    p.rg_solver()
    p.rg_ipb()
    #############################################
    fig = plt.figure()
    # plot
    plt.contour(p.x, p.y, p.z.reshape(p.y.size, p.x.size), 30)
    p.ld.plot(lw = 1.2, ms = 1)
    p.so.plot(lw = 1.2, ms = 1)
    plt.axis('equal')
    plt.axis('square')
    plt.show(block=False)
    if savefig:
      print('sii')
      plt.savefig('runs/fig-thesis/%s_rg_%s%s_alpha%s_no%snd%snb%s_radius%s.eps' %(prename, name, k, p.alpha, no, nd, p.sb.n, 6), bbox_inches='tight')
    return

#######################
def test_init(p=()):
  if p == ():
    p = m.EIT()
  else:
    print("imported p")
  p.domain('one_ellipse')
  p.meshgrid()
  p.theta = 0
  return p;

def test(p=()):
  if p == ():
    p = m.EIT()
  else:
    print("imported p")
  p.domain()
  p.rg_solver()
  p.rg_ipb()
  p.plot()
  inp = input('Press')
  # p.solver()
  # p.ipb()
  # p.plot()
  # inp = input('Press')
  return p
def test_ipb(p=()):
  if p == ():
    p = m.EIT()
  else:
    print("imported p")
  p.domain()
  p.sb = sg.Segment(80, f_inargs = (sh.circle, (0, 6)), quad='ps')
  p.alpha = 1e-15
  p.meshgrid((-2, 2, 60))
  p.theta = 0
  p.rg_solver()
  p.rg_ipb()
  p.plot()
  p.solver()
  p.ipb()
  p.plot()
  inp = input('Press')
  return p

  
if __name__ == "__main__":
  # first_compare()
  # first_compare(gm.sg_one_triangle(nd), 'first_compare', 'triangle')
  # first_compare(gm.one_kite(nd), 'first_compare', 'kite')
  # first_compare(gm.sg_one_drop(nd), 'first_compare', 'drop')
  test_ipb()
  end = input('Press enter')
