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
  f_ex2d, s2d = solve_inc(800, gms)
  for n in rng:
    f2d, s = solve_inc(n, gms)
    f_ex_n = v_ex(s.x)
    ###########################################
    new_err = max(abs(f2d - f_ex2d))
    ############################################
    err = np.concatenate((err, [new_err]))
  return err

def solve_inc(n, gms):
  s = gms(n)
  v_p = ly.scalar(v_p_ex(s.x), s.nx)
  so = s
  #############
  # compute L0
  ld = gm.sg_one_kite(int(3/4 * s.n))
  allpsi = dpb.mapNtoD(ld = ld, lo = so, g = v_p, c = 3)
  Lo = ly.layerpotS(s=so)
  Ld = ly.layerpotS(s=ld, t=so)
  L = Lo.dot(allpsi[0:so.n]) + Ld.dot(allpsi[so.n::])
  onemean = np.ones(so.n).dot(np.diagflat(so.w).dot(L)) / sum(so.w)
  # L0 = L0 - np.array([means for k in range(so.n)])
  #############
  mp = refined.normErrInf()
  Lo = ly.layerpotS(s=so, t=mp)
  Ld = ly.layerpotS(s=ld, t=mp)
  v = Lo.dot(allpsi[0:so.n]) + Ld.dot(allpsi[so.n::])
  v = v - onemean
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

def thesis(n_ex=199):
  rng, err = total(int(n_ex / 2) * 2, gm.sg_one_ellipse)
  # fig = test_plots.plot_loglogscale(rng, err)
  fig = plt.figure()
  plt.plot(rng, err, 'k+-', lw=1, ms=4, ls=':')
  ax = fig.add_subplot(111)
  ax.set_yscale('log')
  # ax.set_xscale('log')
  pnt = ((rng[-1]) * 0.6, (err[-1]))
  ax.annotate('error = %s' % np.float32(err[-1]), xy=pnt , textcoords='data')
  plt.xlabel('n')
  plt.ylabel('log(err)')
  plt.title('Ellipse')
  # plt.axis('square')
  plt.show(block=False)
  if savefig:
    plt.savefig('runs/fig-thesis/convergence_inclusion_one_ellipse_one_kite_inf.eps', bbox_inches='tight')
  return


if __name__ == '__main__':
  # n_ex = sys.argv[1]
  thesis(int(202))
  ret = input('Press')
