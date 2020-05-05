from src import *
import tools_plots

savefig = False

def v_ex(z):
  return z.real**3 - 3 * z.real * z.imag**2
def v_p_ex(z):
  return (3 * z.real**2 - 3 * z.imag**2) + 1j * ( - 6 * z.real * z.imag)
###########################################################################################################
def sg_so(n):
  return sg.Segment(n, f_inargs=(sh.circle, (0, 3)), quad='ps')
###########################################################################################################
def solve(n, gms):
  s = gms(n)
  v_p = ly.scalar(v_p_ex(s.x), s.nx)
  v = ipb.computeL0(so = s, T = v_p)
  return v, s
def solve_inc(n, gms):
  s = gms(n)
  v_p = ly.scalar(v_p_ex(s.x), s.nx)
  v = ipb.computeL(ld = gm.sg_one_kite(int(2/3 * s.n)), so = s, T = v_p, c=3)
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
  v = computeL_exact(ld = (), so = s, T = (), c=3, allpsi=allpsi)
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
    # new_err = max(abs(f - f_ex_n))
    

    f_fem = np.loadtxt('runs/ff++/lapu-%s.txt' %n)
    f_fem = f_fem - sum(f_fem * s.w) / sum(s.w)
    new_fem_err = np.sqrt(sum(s.w * (f_fem - f_ex_n)**2))
    # new_fem_err = max(abs(f_fem - f_ex_n))

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
  for n in rng:
    f, s = solve_inc(n, gms)
    f_ex_n = np.loadtxt('runs/ff++/un-%s.txt' %n)    
    # f_ex_n, s2 = solve_inc_exact(n, gms, allpsi)    
    f_ex_n = f_ex_n - sum(f_ex_n * s.w) / sum(s.w)

    new_err = np.sqrt(sum(s.w * (f - f_ex_n)**2))
    # new_err = max(abs(f - f_ex_n))
    
    f_fem = np.loadtxt('runs/ff++/incu-%s.txt' %n)
    f_fem = f_fem - sum(f_fem * s.w) / sum(s.w)
    new_fem_err = np.sqrt(sum(s.w * (f_fem - f_ex_n)**2))
    # new_fem_err = max(abs(f_fem - f_ex_n))

    err = np.concatenate((err, [new_err]))
    fem_err = np.concatenate((fem_err, [new_fem_err]))
    # plt.plot(s_ex.t, f_ex, '+-')
    # plt.plot(s.t, f_fem, '+-', color = '%s'% float(n/rng[-1]), ms = '0.3')
    # plt.plot(s_ex.t, A.dot(f),'+-')
  # plt.show(block=False)
  # end = input('Press')
  return err, fem_err

def thesis(n_ex = 202):
  rng = range(10 + (n_ex % 2), n_ex, 10)
  err, fem_err = iters(rng, sg_so)
  fig = tools_plots.plot_loglogscale(rng, err)
  ax = fig.add_subplot(111)
  pnt = ((rng[-1]), (err[-1]))
  # plt.plot(pnt[0], pnt[1],'kp')
  ax.annotate('error = %s' % np.float32(err[-1]), xy=pnt , textcoords='data')
  plt.xlabel('log(n)')
  plt.ylabel('log(err)')
  plt.show(block=False)

  fig = tools_plots.plot_loglogscale(rng, fem_err)
  ax = fig.add_subplot(111)
  pnt = ((rng[-1]), (fem_err[-1]))
  # plt.plot(pnt[0], pnt[1],'kp')
  ax.annotate('error = %s' % np.float32(fem_err[-1]), xy=pnt , textcoords='data')
  plt.xlabel('log(n)')
  plt.ylabel('log(fem err)')
  plt.show(block=False)
  #################################################################################
  err, fem_err = iters_inc(rng, sg_so)
  fig = tools_plots.plot_loglogscale(rng, err)
  ax = fig.add_subplot(111)
  pnt = ((rng[-1] ** 0.6), (err[-1]))
  # plt.plot(pnt[0], pnt[1],'kp')
  ax.annotate('error = %s' % np.float32(err[-1]), xy=pnt , textcoords='data')
  plt.xlabel('log(n)')
  plt.ylabel('log(err)')
  plt.show(block=False)
  if savefig:
    plt.savefig('runs/fig-thesis/convergence_inclusion_onekite_bemfem_bem.eps', bbox_inches='tight')
  #############
  fig = tools_plots.plot_loglogscale(rng, fem_err)
  ax = fig.add_subplot(111)
  pnt = ((rng[-1] ** 0.8), (fem_err[-1]))
  # plt.plot(pnt[0], pnt[1],'kp')
  ax.annotate('error = %s' % np.float32(fem_err[-1]), xy=pnt , textcoords='data')
  plt.xlabel('log(n)')
  plt.ylabel('log(fem err)')
  plt.show(block=False)
  if savefig:
    plt.savefig('runs/fig-thesis/convergence_inclusion_onekite_bemfem_fem.eps', bbox_inches='tight')
  return
if __name__ == "__main__":
  thesis()
  ret = input("Press")
