from src import *
import test_plots

savefig = False

ellipses = (
    sg.Layer([sg.Boundary([sg.Segment(80, f_inargs = (sh.ellipse, (0, 2, 1)), quad='ps') ])]), 
    sg.Layer([sg.Boundary([sg.Segment(80, f_inargs = (sh.ellipse, (0, 1.8, 1.2)), quad='ps', aff=(1, np.exp(1j * np.pi/4))) ])])
)

def first_compare(ld = (), prename = (), name = () ):
  if ld == ():
    n = len(ellipses)
    prename, name = 'one_ellipse', 'ellipse'
  else:
    n = 1
  for k in range(0, n):
    p = m.EIT()
    p.domain('one_ellipse', nso = 150)
    # p.ld = sg.Layer([sg.Boundary([sg.Segment(40, f_inargs = (sh.ellipse, (0, 2, 1)), quad='ps') ])])
    # p.ld = sg.Layer([sg.Boundary([sg.Segment(80, f_inargs = (sh.ellipse, (0, 1.5, 1)), quad='ps', aff=(1, 1.2*np.exp(1j * np.pi/4))) ])])
    if ld == ():
      p.ld = ellipses[k]
    else:
      p.ld = ld
    p.meshgrid((-3, 3, 120))
    p.m0 = 40
    p.noiselevel = 0
    #############################################
    # FACTORIZATION METHOD
    ##############################################
    p.fact_solver()
    p.fact_ieig()
    #############################################
    fig = plt.figure()
    # plot
    z = np.array(p.z).reshape((p.x.size, p.y.size))
    fig = plt.contourf(p.x, p.y, z, 1)
    path = fig.collections[1].get_paths()[0]
    v = path.vertices
    # plt.plot(v[:,0], v[:,1],'r-', lw=1.5)
    p.ld.plot(lw = 1.5, ms = 1)
    plt.axis('equal')
    plt.axis('square')
    plt.show(block=False)
    if savefig:
      plt.savefig('runs/fig-thesis/%s_fm_%s%s.eps' %(prename, name, k), bbox_inches='tight')
  return
if __name__ == "__main__":
  # first_compare()
  # first_compare(gm.sg_one_triangle(80), 'first_compare', 'triangle')
  first_compare(gm.sg_one_kite(80), 'first_compare', 'kite')
  # first_compare(gm.sg_one_drop(80), 'first_compare', 'drop')
  end = input('Press enter')
