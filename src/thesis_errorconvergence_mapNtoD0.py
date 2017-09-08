from src import *
import test_plots

savefig = False

def v_ex(z):
  return z.real**3 - 3 * z.real * z.imag**2
def v_p_ex(z):
  return (3 * z.real**2 - 3 * z.imag**2) + 1j * ( - 6 * z.real * z.imag)
###########################################################################################################
def e0(tt, st0, st1):
  # first element function
  return (tt - st0) / (st1 - st0)
def e1(tt, st0, st1):
  # second element function
  return (st1 - tt) / (st1 - st0)
def truncate(f):
  return f * (f>=0) * (f<=1)
def truncate_up(f):
  # truncate in the exact quadrature nodes the sum of two element functions
  return 1 + (f - 1) * (f<=1)
def elem_deg1(s, t):
  # compute all element functions for sources s and targets t
  st_ext = np.concatenate( ([s.t[-1] - 1], s.t, [s.t[0] + 1]) )
  A = np.empty((len(t.t), len(s.t)))
  for j in range(1, len(st_ext) - 1 ):
    A[:, j-1] = truncate_up(truncate(e0(t.t, st_ext[j-1], st_ext[j])) + truncate(e1(t.t, st_ext[j], st_ext[j+1])))\
              + truncate_up(truncate(e0(t.t - 1, st_ext[j-1], st_ext[j])) + truncate(e1(t.t - 1, st_ext[j], st_ext[j+1])))\
              + truncate_up(truncate(e0(t.t + 1, st_ext[j-1], st_ext[j])) + truncate(e1(t.t + 1, st_ext[j], st_ext[j+1])))
  return A  
###########################################################################################################
def iters(rng, f_ex, s_ex, gms):
  err = []
  for n in rng:
    f, s = solve(n, gms)
    # A = elem_deg1(s, s_ex)
    # new_err = linalg.norm((A.dot(f) - f_ex) * s_ex.w)
    f_ex_n = v_ex(s.x)
    f_ex_n = f_ex_n - sum(f_ex_n * s.w) / sum(s.w)
    new_err = np.sqrt(sum(s.w * (f - f_ex_n)**2))
    # new_err = max(abs(f - f_ex_n))
    err = np.concatenate((err, [new_err]))
    # plt.plot(s_ex.t, f_ex, '+-')
    # plt.plot(s.t, f, '+-')
    # plt.plot(s_ex.t, A.dot(f),'+-')
    # plt.show(block=False)
    # end = input('Press')
  return err

def solve(n, gms):
  # s = gm.sg_one_kite(n)
  s = gms(n)
  v_p = ly.scalar(v_p_ex(s.x), s.nx)
  v = ipb.computeL0(so = s, T = v_p)
  # L0 = ipb.computeL0(so = s, T = np.eye(s.n))
  # return L0.dot(v_p)
  return v, s

def exact(n, gms):
  # s_ex = gm.sg_one_kite(n)
  s_ex = gms(n)
  f_ex = v_ex(s_ex.x)
  f_ex = f_ex - sum(f_ex * s_ex.w) / sum(s_ex.w)
  return f_ex, s_ex
  

def total(n_ex, gms=gm.sg_one_triangle):
  f_ex, s_ex = exact(n_ex, gms)
  rng = range(10 + (n_ex % 2), n_ex, 10)
  err = iters(rng, f_ex, s_ex, gms)
  # plt.plot(err, '+-')
  # plt.show(block=False)
  return rng, err
def thesis(n_ex=199):
  rng, err = total(int(n_ex / 2) * 2 + 1, gm.sg_one_triangle)
  fig = test_plots.plot_loglogscale(rng, err)
  ax = fig.add_subplot(111)
  pnt = ((rng[-1]), (err[-1]))
  # plt.plot(pnt[0], pnt[1],'kp')
  ax.annotate('error = %s' % np.float32(err[-1]), xy=pnt , textcoords='data')
  plt.xlabel('log(n)')
  plt.ylabel('log(err)')
  plt.title('Triangle')
  plt.show(block=False)
  if savefig:
    plt.savefig('runs/fig-thesis/convergence_laplace_one_triangle.eps', bbox_inches='tight')
  ##############################
  rng, err = total(int(n_ex / 2) * 2, gm.sg_one_ellipse)
  fig = test_plots.plot_loglogscale(rng, err)
  ax = fig.add_subplot(111)
  pnt = ((rng[-1]), (err[-1]))
  # plt.plot(pnt[0], pnt[1],'kp')
  ax.annotate('error = %s' % np.float32(err[-1]), xy=pnt , textcoords='data')
  plt.xlabel('log(n)')
  plt.ylabel('log(err)')
  plt.title('Ellipse')
  plt.show(block=False)
  if savefig:
    plt.savefig('runs/fig-thesis/convergence_laplace_one_ellipse.eps', bbox_inches='tight')
  ##############################
  rng, err = total(int(n_ex / 2) * 2, gm.sg_one_kite)
  fig = test_plots.plot_loglogscale(rng, err)
  ax = fig.add_subplot(111)
  pnt = ((rng[-1]), (err[-1]))
  # plt.plot(pnt[0], pnt[1],'kp')
  ax.annotate('error = %s' % np.float32(err[-1]), xy=pnt , textcoords='data')
  plt.xlabel('log(n)')
  plt.ylabel('log(err)')
  plt.title('Kite')
  plt.show(block=False)
  if savefig:
    plt.savefig('runs/fig-thesis/convergence_laplace_one_kite.eps', bbox_inches='tight')
  #############################
  rng, err = total(int(n_ex / 2) * 2 + 1, gm.sg_one_drop)
  fig = test_plots.plot_loglogscale(rng, err)
  ax = fig.add_subplot(111)
  pnt = ((rng[-1]), (err[-1]))
  # plt.plot(pnt[0], pnt[1],'kp')
  ax.annotate('error = %s' % np.float32(err[-1]), xy=pnt , textcoords='data')
  plt.xlabel('log(n)')
  plt.ylabel('log(err)')
  plt.title('Drop')
  plt.show(block=False)
  if savefig:
    plt.savefig('runs/fig-thesis/convergence_laplace_one_drop.eps', bbox_inches='tight')
  return

if __name__ == '__main__':
  n_ex = sys.argv[1]
  thesis(int(n_ex))
  ret = input('Press')
