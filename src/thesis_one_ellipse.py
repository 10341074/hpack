from src import *
import test_plots
# import dill

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
    p.domain('one_ellipse', nso = 100)
    # p.ld = sg.Layer([sg.Boundary([sg.Segment(40, f_inargs = (sh.ellipse, (0, 2, 1)), quad='ps') ])])
    # p.ld = sg.Layer([sg.Boundary([sg.Segment(80, f_inargs = (sh.ellipse, (0, 1.5, 1)), quad='ps', aff=(1, 1.2*np.exp(1j * np.pi/4))) ])])
    if ld == ():
      p.ld = ellipses[k]
    else:
      p.ld = ld
    p.meshgrid((-3, 3, 80))
    p.noiselevel = 0
    p.solver()
    p.ipb()
    # p.plot()
    p.plot_pre()
    fig = plt.figure()
    plt.contour(p.x, p.y, p.z.reshape(p.y.size, p.x.size), 30)
    p.ld.plot(lw = 0.8, ms = 1)
    p.so.plot(lw = 0.8, ms = 1)
    plt.axis('equal')
    plt.axis('square')
    # plt.show(block=False)
    plt.savefig('runs/fig-thesis/%s_lsm_%s%s.eps' %(prename, name, k), bbox_inches='tight')
    # test_plots.plot_contourf(p, ())

    # alpha = np.array([1e-16 + 1e-16 * k for k in range(10)])
    # alpha = np.concatenate((alpha, [1e-15*10**k for k in range(8)]))
    #p.test_alpha(alpha=alpha)

    # L0 = computeL0(p.so, p.so.BX, p.pp)
    # z = L0.dot(p.isolver.save_sol[0])
    # p.plot(z)
    ###########################
    # p.solver()
    # it_alpha = 8
    # p.ipb_opt(it_alpha=it_alpha)
    # p.alpha_fixed_ratio(0)
    # p.plot()
    # p.alpha_fixed_ratio(-1)
    # p.plot()
    #############################
    # dill.dump_session('dill.pkl')
    # dill.load_session('dill.pkl')
    # CORRECT
    ###########################
    plt.figure()
    p.fact_solver()
    p.fact_ieig()
    v = test_plots.plot_contourf_1(p, ())
    # p.plot()
    plt.plot(v[:,0], v[:,1],'r-')
    o = test_plots.stats(v)
    plt.plot([o[0].real, v[o[3],0]], [o[0].imag, v[o[3],1]], 'b-')
    print(o)
    p.ld.plot(lw = 0.8, ms = 1)
    plt.axis('equal')
    plt.axis('square')
    # plt.show(block=False)
    plt.savefig('runs/fig-thesis/%s_fm_%s%s.eps' %(prename, name, k), bbox_inches='tight')
  return
if __name__ == "__main__":
  first_compare()
  # first_compare(gm.sg_one_triangle(80), 'first_compare', 'triangle')
  # first_compare(gm.sg_one_kite(80), 'first_compare', 'kite')
  # first_compare(gm.sg_one_drop(80), 'first_compare', 'drop')
  end = input('Press enter')

