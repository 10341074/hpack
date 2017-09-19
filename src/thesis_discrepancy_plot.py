from src import *
import test_plots
import dill

savefig = False

def plot_compare(filename = (), prename = (), name = () ):
  n = 1
  for k in range(0, n):
    dill.load_session(filename)
    for it_alpha_fixed_ratio in np.array([0, 14, 15, 19]):
      p.alpha_fixed_ratio(it_alpha_fixed_ratio)
      # p.plot()
      p.plot_pre()
      fig = plt.figure()
      plt.contour(p.x, p.y, p.z.reshape(p.y.size, p.x.size), 30)
      p.ld.plot(lw = 0.8, ms = 1)
      p.so.plot(lw = 0.8, ms = 1)
      plt.colorbar()
      plt.axis('equal')
      plt.axis('square')
      plt.show(block=False)
      if savefig:
        plt.savefig('runs/fig-thesis/%s_lsm_%s%s_it_alpha_%s.eps' %(prename, name, k, it_alpha_fixed_ratio), bbox_inches='tight')
    ##########################################################################################################################
    it_alpha_fixed_ratio = 19
    p.z = p.isolver.save_alpha[:, it_alpha_fixed_ratio]
    p.plot_pre()
    fig = plt.figure()
    plt.contour(p.x, p.y, p.z.reshape(p.y.size, p.x.size), 30)
    p.ld.plot(lw = 0.8, ms = 1)
    p.so.plot(lw = 0.8, ms = 1)
    plt.colorbar()
    plt.axis('equal')
    plt.axis('square')
    plt.title('alpha')
    plt.show(block=False)
    if savefig:
      plt.savefig('runs/fig-thesis/%s_lsm_%s%s_alpha_it_alpha_%s.eps' %(prename, name, k, it_alpha_fixed_ratio), bbox_inches='tight')
    # ###########################################################################################################################
    # # FACTORIZATION
    # ##########################################################################################################################
    # plt.figure()
    # p.fact_solver()
    # p.fact_ieig()
    # v = test_plots.plot_contourf_1(p, ())
    # # p.plot()
    # plt.plot(v[:,0], v[:,1],'r-')
    # o = test_plots.stats(v)
    # plt.plot([o[0].real, v[o[3],0]], [o[0].imag, v[o[3],1]], 'b-')
    # print(o)
    # p.ld.plot(lw = 0.8, ms = 1)
    # plt.axis('equal')
    # plt.axis('square')
    # plt.show(block=False)
    # if savefig:
    #   plt.savefig('runs/fig-thesis/%s_fm_%s%s.eps' %(prename, name, k), bbox_inches='tight')
  return

def plot_alpha(filename = (), prename = (), name = () ):
  dill.load_session(filename)
  k = 0
  it_alpha_fixed_ratio = 0
  p.z = p.isolver.save_alpha[:, it_alpha_fixed_ratio]
  p.plot_pre()
  fig = plt.figure()
  plt.contour(p.x, p.y, p.z.reshape(p.y.size, p.x.size), 30)
  p.ld.plot(lw = 0.8, ms = 1)
  p.so.plot(lw = 0.8, ms = 1)
  plt.colorbar()
  plt.title('alpha')
  plt.axis('equal')
  plt.axis('square')
  plt.show(block=False)
  if savefig:
    plt.savefig('runs/fig-thesis/%s_lsm_%s%s_alpha_it_alpha_%s.eps' %(prename, name, k, it_alpha_fixed_ratio), bbox_inches='tight')

if __name__ == "__main__":
  # plot_compare('runs/dill/dill_two_ellipse_delta1e-6.pkl', 'discrepancy', 'two_ellipse_delta1e-6')
  # plot_compare('runs/dill/dill_three_ellipse_delta1e-6.pkl', 'discrepancy', 'three_ellipse_delta1e-6')
  # plot_compare('runs/dill/dill_sharp_ellipse_delta1e-6.pkl', 'discrepancy', 'sharp_ellipse_delta1e-6')
  # plot_compare('runs/dill/dill_two_ellipse_delta1e-6_app.pkl', 'discrepancy', 'two_ellipse_delta1e-6')
  # plot_compare('runs/dill/dill_three_ellipse_delta1e-6_app.pkl', 'discrepancy', 'three_ellipse_delta1e-6')
  # plot_compare('runs/dill/dill_sharp_ellipse_delta1e-6_app.pkl', 'discrepancy', 'sharp_ellipse_delta1e-6')
  plot_compare('runs/dill/dill_three_ellipse_noiseK1e-3_delta1e-6_bis.pkl', 'discrepancy_noiseK1e-3', 'three_ellipse_delta1e-6')
  end = input('Press enter')
