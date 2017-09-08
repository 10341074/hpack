from src import *
import test_plots

savefig = False

def v_ex(z):
  return z.real**3 - 3 * z.real * z.imag**2
def v_p_ex(z):
  return (3 * z.real**2 - 3 * z.imag**2) + 1j * ( - 6 * z.real * z.imag)
###########################################################################################################

def sg_so(n):
  return sg.Segment(n, f_inargs=(sh.circle, (0, 3)), quad='p')
def solve(n, gms):
  s = gms(n)
  v_p = ly.scalar(v_p_ex(s.x), s.nx)
  v = ipb.computeL0(so = s, T = v_p)
  # L0 = ipb.computeL0(so = s, T = np.eye(s.n))
  # return L0.dot(v_p)
  return v, s
def solve_inc(n, gms):
  s = gms(n)
  v_p = ly.scalar(v_p_ex(s.x), s.nx)
  v = ipb.computeL(ld = gm.sg_one_kite(int(3/4 * s.n)), so = s, T = v_p, c=3)
  # L0 = ipb.computeL0(so = s, T = np.eye(s.n))
  # return L0.dot(v_p)
  return v, s
def computeL_allpsi_exact(ld, so, T, c):
  print('computing L EXACT')
  so_ex = sg_so(500)
  T = ly.scalar(v_p_ex(so_ex.x), so_ex.nx)
  allpsi = dpb.mapNtoD(so_ex, ld, T, c, so_ex.s0)
  # allpsi = dpb.mapNtoDD0(so, ld, T, c, so.s0)
  # allpsi = dpb.mapNtoDD_correctedinfirst(so, ld, T, c, so.s0)
  # allpsi = dpb.mapNtoDD_left(so, ld, T, c, so.s0)
  return allpsi
def computeL_exact(ld, so, T, c, allpsi):  
  so_ex = sg_so(500)
  ld_ex = gm.sg_one_kite(int(2/3 * 500))
  Lo = ly.layerpotS(s=so_ex, t=so)
  Ld = ly.layerpotS(s=ld_ex, t=so)
  L = Lo.dot(allpsi[0:so_ex.n]) + Ld.dot(allpsi[so_ex.n::])
  # means = sum(np.diagflat(so.w).dot(L)) / sum(so.w) # correct? strange sum by rows
  means = np.ones(so.n).dot(np.diagflat(so.w).dot(L)) / sum(so.w)
  L = L - np.array([means for k in range(so.n)])
  return L

def solve_inc_exact(n, gms, allpsi):
  s = gms(n)
  # v_p = ly.scalar(v_p_ex(s.x), s.nx)
  # v = computeL_exact(ld = gm.sg_one_kite(int(2/3 * s.n)), so = s, T = v_p, c=3, allpsi=allpsi)
  v = computeL_exact(ld = (), so = s, T = (), c=3, allpsi=allpsi)
  # L0 = ipb.computeL0(so = s, T = np.eye(s.n))
  # return L0.dot(v_p)
  return v, s

def iters(rng, gms):
  err, fem_err = [], []
  for n in rng:
    f, s = solve(n, gms)
    # A = elem_deg1(s, s_ex)
    # new_err = linalg.norm((A.dot(f) - f_ex) * s_ex.w)
    f_ex_n = v_ex(s.x)
    f_ex_n = f_ex_n - sum(f_ex_n * s.w) / sum(s.w)
    new_err = np.sqrt(sum(s.w * (f - f_ex_n)**2))
    
    f_fem = np.loadtxt('runs/ff++/lapu-%s.txt' %n)
    f_fem = f_fem - sum(f_fem * s.w) / sum(s.w)
    new_fem_err = np.sqrt(sum(s.w * (f_fem - f_ex_n)**2))

    err = np.concatenate((err, [new_err]))
    fem_err = np.concatenate((fem_err, [new_fem_err]))
    # plt.plot(s_ex.t, f_ex, '+-')
    # plt.plot(s.t, f_fem, '+-', color = '%s'% float(n/rng[-1]), ms = '0.3')
    # plt.plot(s_ex.t, A.dot(f),'+-')
  # plt.show(block=False)
  # end = input('Press')
  return err, fem_err
def iters_inc(rng, gms):
  err, fem_err = [], []
  # allpsi = computeL_allpsi_exact(ld = gm.sg_one_kite(int(2/3*500)), so = (), T = (), c=3)
  f_ex_n, s2 = solve_inc(800, gms)
  print("FIRST = ", s2.x[0])
  f_ex_n = f_ex_n - sum(f_ex_n * s2.w) / sum(s2.w)
  for n in rng:
    f, s = solve_inc(n, gms)
    print("FIRST = ", s.x[0])
    ind = np.arange(s.n) * 800 / s.n
    f_extr = np.array([f_ex_n[int(k)] for k in ind])
    # A = elem_deg1(s, s_ex)
    # new_err = linalg.norm((A.dot(f) - f_ex) * s_ex.w)

    # f_ex_n = np.loadtxt('runs/ff++/un-%s.txt' %n)    
    # f_ex_n, s2 = solve_inc_exact(n, gms, allpsi)
    
    # f_ex_n = f_ex_n - sum(f_ex_n * s.w) / sum(s.w)

    new_err = np.sqrt(sum(s.w * (f - f_extr)**2))
    err = np.concatenate((err, [new_err]))
    
    # plt.plot(s_ex.t, f_ex, '+-')
    # plt.plot(s.t, f_fem, '+-', color = '%s'% float(n/rng[-1]), ms = '0.3')
    # plt.plot(s_ex.t, A.dot(f),'+-')
  # plt.show(block=False)
  # end = input('Press')
  return err

# def exact(n, gms):
#   s_ex = gms(n)
#   f_ex = v_ex(s_ex.x)
#   f_ex = f_ex - sum(f_ex * s_ex.w) / sum(s_ex.w)
#   return f_ex, s_ex
  

# def total(n_ex, gms=gm.sg_one_triangle):
#   # f_ex, s_ex = exact(n_ex, gms)
#   rng = range(10 + (n_ex % 2), n_ex, 10)
#   err, fem_err = iters_inc(rng, gms)
#   return rng, err, fem_err

def thesis(n_ex = 202):
  # rng = range(10 + (n_ex % 2), n_ex, 10)
  # err, fem_err = iters(rng, sg_so)
  # fig = test_plots.plot_loglogscale(rng, err)
  # ax = fig.add_subplot(111)
  # pnt = ((rng[-1]), (err[-1]))
  # # plt.plot(pnt[0], pnt[1],'kp')
  # ax.annotate('error = %s' % np.float32(err[-1]), xy=pnt , textcoords='data')
  # plt.xlabel('log(n)')
  # plt.ylabel('log(err)')
  # plt.show(block=False)

  # fig = test_plots.plot_loglogscale(rng, fem_err)
  # ax = fig.add_subplot(111)
  # pnt = ((rng[-1]), (fem_err[-1]))
  # # plt.plot(pnt[0], pnt[1],'kp')
  # ax.annotate('error = %s' % np.float32(fem_err[-1]), xy=pnt , textcoords='data')
  # plt.xlabel('log(n)')
  # plt.ylabel('log(fem err)')
  # plt.show(block=False)
  #################################################################################
  rng = np.array([10, 16, 20, 32, 40, 50, 80, 100, 160, 200, 400])
  err = iters_inc(rng, sg_so)
  fig = test_plots.plot_loglogscale(rng, err)
  ax = fig.add_subplot(111)
  pnt = ((rng[-1]), (err[-1]))
  # plt.plot(pnt[0], pnt[1],'kp')
  ax.annotate('error = %s' % np.float32(err[-1]), xy=pnt , textcoords='data')
  plt.xlabel('log(n)')
  plt.ylabel('log(err)')
  plt.show(block=False)
  if savefig:
    plt.savefig('runs/fig-thesis/convergence_inclusion_onekite_bemfem_bem_bemerr.eps', bbox_inches='tight')
  return
if __name__ == "__main__":
  thesis()
  ret = input("Press")
