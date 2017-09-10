from src import *
import test_plots
import refined

savefig = False

def v_ex(z):
  return z.real**3 - 3 * z.real * z.imag**2
def v_p_ex(z):
  return (3 * z.real**2 - 3 * z.imag**2) + 1j * ( - 6 * z.real * z.imag)
###########################################################################################################
def iters(rng, f_ex, s_ex, gms):
  err = []
  sbig = gms(600)
  for n in rng:
    f, s = solve(n, gms)
    f_ex_n = v_ex(s.x)
    f_ex_n = f_ex_n - sum(f_ex_n * s.w) / sum(s.w)
    # new_err = np.sqrt(sum(s.w * (f - f_ex_n)**2))
    ############################################
    # error with P1 definition on fitted mesh
    # new_err = refined.normErr(gms, (), (), v_ex, f, sbig, s, 'Inf')
    ###########################################
    new_err = max(abs(f - f_ex_n))
    err = np.concatenate((err, [new_err]))
  return err

def solve(n, gms):
  s = gms(n)
  v_p = ly.scalar(v_p_ex(s.x), s.nx)
  v = ipb.computeL0(so = s, T = v_p)
  return v, s

def exact(n, gms):
  s_ex = gms(n)
  f_ex = v_ex(s_ex.x)
  f_ex = f_ex - sum(f_ex * s_ex.w) / sum(s_ex.w)
  return f_ex, s_ex

def total(n_ex, gms=gm.sg_one_triangle):
  f_ex, s_ex = exact(n_ex, gms)
  rng = range(10 + (n_ex % 2), n_ex, 10)
  err = iters(rng, f_ex, s_ex, gms)
  return rng, err

def thesis(n_ex=202):
  #############################
  # rng, err = total(int(n_ex / 2) * 2 + 1, gm.sg_one_triangle)
  # fig = test_plots.plot_loglogscale(rng, err)
  # ax = fig.add_subplot(111)
  # pnt = ((rng[-1] ** 0.8), (err[-1]))
  # # plt.plot(pnt[0], pnt[1],'kp')
  # ax.annotate('error = %s' % np.float32(err[-1]), xy=pnt , textcoords='data')
  # plt.xlabel('log(n)')
  # plt.ylabel('log(err)')
  # plt.title('Triangle')
  # plt.show(block=False)
  # if savefig:
  #   plt.savefig('runs/fig-thesis/convergence_laplace_one_triangle.eps', bbox_inches='tight')
  ##############################
  rng, err = total(int(n_ex / 2) * 2, gm.sg_one_ellipse)
  fig = test_plots.plot_loglogscale(rng, err)
  ax = fig.add_subplot(111)
  pnt = ((rng[-1] ** 0.8), (err[-1] ** 0.9))
  # plt.plot(pnt[0], pnt[1],'kp')
  ax.annotate('error = %s' % np.float32(err[-1]), xy=pnt , textcoords='data')
  plt.xlabel('log(n)')
  plt.ylabel('log(err)')
  plt.title('Ellipse')
  plt.show(block=False)
  if savefig:
    plt.savefig('runs/fig-thesis/convergence_laplace_one_ellipse_pointwise.eps', bbox_inches='tight')
  # ###############################
  # rng, err = total(int(n_ex / 2) * 2, gm.sg_one_kite)
  # fig = test_plots.plot_loglogscale(rng, err)
  # ax = fig.add_subplot(111)
  # pnt = ((rng[-1] ** 0.8), (err[-1]))
  # # plt.plot(pnt[0], pnt[1],'kp')
  # ax.annotate('error = %s' % np.float32(err[-1]), xy=pnt , textcoords='data')
  # plt.xlabel('log(n)')
  # plt.ylabel('log(err)')
  # plt.title('Kite')
  # plt.show(block=False)
  # if savefig:
  #   plt.savefig('runs/fig-thesis/convergence_laplace_one_kite.eps', bbox_inches='tight')
  # #############################
  # rng, err = total(int(n_ex / 2) * 2 + 1, gm.sg_one_drop)
  # fig = test_plots.plot_loglogscale(rng, err)
  # ax = fig.add_subplot(111)
  # pnt = ((rng[-1] ** 0.8), (err[-1]))
  # # plt.plot(pnt[0], pnt[1],'kp')
  # ax.annotate('error = %s' % np.float32(err[-1]), xy=pnt , textcoords='data')
  # plt.xlabel('log(n)')
  # plt.ylabel('log(err)')
  # plt.title('Drop')
  # plt.show(block=False)
  # if savefig:
  #   plt.savefig('runs/fig-thesis/convergence_laplace_one_drop.eps', bbox_inches='tight')
  return


if __name__ == '__main__':
  # n_ex = sys.argv[1]
  thesis(int(202))
  ret = input('Press')
