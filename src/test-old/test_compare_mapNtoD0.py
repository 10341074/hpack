from src import *
import tools_plots

def v_ex(z):
  return z.real**3 - 3 * z.real * z.imag**2
def v_p_ex(z):
  return (3 * z.real**2 - 3 * z.imag**2) + 1j * ( - 6 * z.real * z.imag)
#######################################################################
def iters(rng, f_ex, s_ex, gms):
  err = []
  for n in rng:
    f, s = solve(n, gms)
    # A = elem_deg1(s, s_ex)
    # new_err = linalg.norm((A.dot(f) - f_ex) * s_ex.w)
    f_ex_n = v_ex(s.x)
    f_ex_n = f_ex_n - sum(f_ex_n * s.w) / sum(s.w)
    new_err = linalg.norm((f - f_ex_n) * s.w)
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
  # L0 = ipb.computeL0(so = s, T = np.eye(s.n))
  v_p = ly.scalar(v_p_ex(s.x), s.nx)
  v = ipb.computeL0(so = s, T = v_p)
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
  plt.plot(err, '+-')
  plt.show(block=False)
  return rng, err
##################################################################################################
def test_injectivity_N0(so):
  return ipb.computeL0(so = so, T = np.zeros(so.n))
def test_injectivity_ND(so):
  return ipb.computeL(so = so, T = np.zeros(so.n))
def test_injectivity(so, N0=True, ND=False):
  if N0: out(test_injectivity_N0(so), s = 'N0')
  if ND: out(test_injectivity_ND(so), s = 'ND')
  return
def out(v, s):
  e = max(abs(v))
  toll = 1e-10
  out = 'OK' if e < toll else 'KO'
  print(out, ': test injectivity of ', s, ' with max abs e = ', e, ' < ', toll)
  return


def submean(v, w, strname, sub=1):
  print(strname, 'old mean = ', sum(v * w))
  if sub:
    v = v - sum(v*w)/sum(w)
    print('Subtr\n')
  return v

def test_N0(so):
  v_p = ly.scalar(v_p_ex(so.x), so.nx)
  f_ex = v_ex(so.x)
  f_ex = submean(f_ex, so.w, 'Exact')
  e = []
  # direct precond
  v__dp = ly.layerpotS(s=so).dot(dpb.mapNtoD0_left(so, v_p, ()))
  v__dp = submean(v__dp, so.w, 'Direct Precon')
  e = np.concatenate(( e, [test_N0_out(f_ex - v__dp)] ))
  # direct
  v__d = ly.layerpotS(s=so).dot(dpb.mapNtoD0_correctedinfirst(so, v_p, ()))
  v__d = submean(v__d, so.w, 'Direct')
  e = np.concatenate(( e, [test_N0_out(f_ex - v__d)] ))
  # clin precond
  v__clp = ly.layerpotS(s=so).dot(dpb.mapNtoD0_left(so, np.eye(so.n), ())).dot(v_p)
  v__clp = submean(v__clp, so.w, 'Composite Linear Precon')
  e = np.concatenate(( e, [test_N0_out(f_ex - v__clp)] ))
  # clin
  v__cl = ly.layerpotS(s=so).dot(dpb.mapNtoD0_correctedinfirst(so, np.eye(so.n), ())).dot(v_p)
  v__cl = submean(v__cl, so.w, 'Composite Linear')
  e = np.concatenate(( e, [test_N0_out(f_ex - v__cl)] ))
  return e, f_ex, v__dp, v__d, v__clp, v__cl

def test_N0_out(v, strname=()):
  e = max(abs(v))
  toll = 1e-10
  out = 'OK' if e < toll else 'KO'
  print(out, ': test error of ', strname, ' \t with max abs e = ', e, ' < ', toll)
  return e

def test_N0_total(n_ex):
  # f_ex, s_ex = exact(n_ex)
  rng = range(10 + (n_ex % 2), n_ex, 10)
  err = [np.zeros(4)]
  for n in rng:
    s = gm.sg_one_kite(n)
    result = test_N0(s)
    err = np.concatenate(( err, [np.array(result[0])] ))
  for k in range(4):
    plt.plot(err[:, k], '+-')
    plt.show(block=False)
    temp = input('Press')
  return rng, err


if __name__ == '__main__':
  n_ex = int(sys.argv[1])
  test_N0_total(n_ex)
