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
def solve_inc(n, gms):
  s = gms(n)
  v_p = ly.scalar(v_p_ex(s.x), s.nx)
  v = ipb.computeL(ld = gm.sg_one_kite(int(3/4 * s.n)), so = s, T = v_p, c=3)
  return v, s
# def computeL_allpsi_exact(ld, so, T, c):
#   print('computing L EXACT')
#   so_ex = sg_so(500)
#   T = ly.scalar(v_p_ex(so_ex.x), so_ex.nx)
#   allpsi = dpb.mapNtoD(so_ex, ld, T, c, so_ex.s0)
#   # allpsi = dpb.mapNtoDD0(so, ld, T, c, so.s0)
#   # allpsi = dpb.mapNtoDD_correctedinfirst(so, ld, T, c, so.s0)
#   # allpsi = dpb.mapNtoDD_left(so, ld, T, c, so.s0)
#   return allpsi
# def computeL_exact(ld, so, T, c, allpsi):  
#   so_ex = sg_so(500)
#   ld_ex = gm.sg_one_kite(int(2/3 * 500))
#   Lo = ly.layerpotS(s=so_ex, t=so)
#   Ld = ly.layerpotS(s=ld_ex, t=so)
#   L = Lo.dot(allpsi[0:so_ex.n]) + Ld.dot(allpsi[so_ex.n::])
#   # means = sum(np.diagflat(so.w).dot(L)) / sum(so.w) # correct? strange sum by rows
#   means = np.ones(so.n).dot(np.diagflat(so.w).dot(L)) / sum(so.w)
#   L = L - np.array([means for k in range(so.n)])
#   return L

# def solve_inc_exact(n, gms, allpsi):
#   s = gms(n)
#   v = computeL_exact(ld = (), so = s, T = (), c=3, allpsi=allpsi)
#   return v, s

def iters_inc(rng, gms):
  err = []
  # allpsi = computeL_allpsi_exact(ld = gm.sg_one_kite(int(2/3*500)), so = (), T = (), c=3)
  f_ex_n, s2 = solve_inc(800, gms)
  f_ex_n = f_ex_n - sum(f_ex_n * s2.w) / sum(s2.w)

  for n in rng:
    f, s = solve_inc(n, gms)
    ind = np.arange(s.n) * 800 / s.n
    f_extr = np.array([f_ex_n[int(k)] for k in ind])

    new_err = np.sqrt(sum(s.w * (f - f_extr)**2))
    err = np.concatenate((err, [new_err]))
    
    # plt.plot(s_ex.t, f_ex, '+-')
    # plt.plot(s.t, f_fem, '+-', color = '%s'% float(n/rng[-1]), ms = '0.3')
    # plt.plot(s_ex.t, A.dot(f),'+-')
  # plt.show(block=False)
  # end = input('Press')
  return err

def thesis(n_ex = 202):
  rng = np.array([16, 20, 32, 40, 50, 80, 100, 160, 200, 400])
  err = iters_inc(rng, sg_so)

  fig = plt.figure()
  plt.plot(rng, err, 'k+-', lw=1, ms=4, ls=':')
  ax = fig.add_subplot(111)
  ax.set_yscale('log')
  ax.set_xscale('log')
  pnt = ((rng[-1]**0.8), (err[-1]))
  ax.annotate('error = %s' % np.float32(err[-1]), xy=pnt , textcoords='data')
  plt.xlabel('log(n)')
  plt.ylabel('log(err)')
  # plt.axis('square')
  plt.show(block=False)
  if savefig:
    plt.savefig('runs/fig-thesis/convergence_inclusion_onekite_bemfem_bem_bemerr.eps', bbox_inches='tight')
  return
if __name__ == "__main__":
  thesis()
  ret = input("Press")
