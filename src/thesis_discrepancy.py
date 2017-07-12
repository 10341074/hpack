from src import *
import test_plots
import dill

def first_compare(filename = (), prename = (), name = () ):
  n = 1
  for k in range(0, n):
    dill.load_session(filename)
    for it_alpha_fixed_ratio in range(18):
      p.alpha_fixed_ratio(it_alpha_fixed_ratio)
      # p.plot()
      p.plot_pre()
      fig = plt.figure()
      plt.contour(p.x, p.y, p.z.reshape(p.y.size, p.x.size), 30)
      p.ld.plot(lw = 0.8, ms = 1)
      p.so.plot(lw = 0.8, ms = 1)
      plt.axis('equal')
      plt.axis('square')
      # plt.show(block=False)
      plt.savefig('runs/fig-thesis/%s_lsm_%s%s_it_alpha_%s.eps' %(prename, name, k, it_alpha_fixed_ratio), bbox_inches='tight')
    #########################################################
    # test_plots.plot_contourf(p, ())
    #########################################################
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
  first_compare('runs/dill/dill_two_ellipse.pkl', 'discrepancy', 'two_ellipse')
  end = input('Press enter')
