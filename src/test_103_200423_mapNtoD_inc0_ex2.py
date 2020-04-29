from src import *

from data import v1 as v_ex
from data import v1_p as v_p_ex

# --------------------------------------------------------------------------------

def mean(v, s):
  """This function compute weighted mean of v on segment s"""
  return sum(v * s.w) / sum(s.w)
def submean(v, w, strname, sub=1):
  if sub:
    v = v - sum(v*w)/sum(w)
  return v
def test_N0_out(v, strname=()):
  e = max(abs(v))
  toll = 1e-10
  out = 'OK' if e < toll else 'KO'
  print(out, ': test error of ', strname, ' \t with max abs e = ', e, ' < ', toll)
  return e

def test_N0(so):
  v_p = ly.scalar(v_p_ex(so.x), so.nx)
  f_ex = v_ex(so.x)
  f_ex_mean = mean(f_ex, so)
  f_ex = f_ex - f_ex_mean
  e = []
  # --------------------------------------------------------------
  # direct precond
  phi = dpb.mapNtoD0_left(so, v_p, ())
  v__dp = ly.layerpotS(s=so).dot(phi)
  v_mean = mean(v__dp, so)
  v__dp = submean(v__dp, so.w, 'Direct Precon')
  e = np.concatenate(( e, [test_N0_out(f_ex - v__dp, 'Direct Precon')] ))
  # ----------------------
  if False:
    # --- intensity plot ---
    x, y, pp = plot.meshgrid((-3, 3, 50))
    yy = ly.layerpotS(s=so, t=sg.Pointset(x=pp)).dot(phi) - v_mean - (v_ex(pp) - f_ex_mean)
    plot.plot(x, y, yy, vmin=-0.2, vmax=0.2)
    so.plot()
    plt.show(block=False)
  # --------------------------------------------------------------
  # direct
  phi = dpb.mapNtoD0_correctedinfirst(so, v_p, ())
  v__d = ly.layerpotS(s=so).dot(phi)
  v__d = submean(v__d, so.w, 'Direct')
  e = np.concatenate(( e, [test_N0_out(f_ex - v__d, 'Direct')] ))
  # --------------------------------------------------------------
  # clin precond
  phi = dpb.mapNtoD0_left(so, np.eye(so.n), ()).dot(v_p)
  v__clp = ly.layerpotS(s=so).dot(phi)
  v__clp = submean(v__clp, so.w, 'Composite Linear Precon')
  e = np.concatenate(( e, [test_N0_out(f_ex - v__clp, 'Composite Linear Precon')] ))
  # --------------------------------------------------------------
  # clin
  phi = dpb.mapNtoD0_correctedinfirst(so, np.eye(so.n), ()).dot(v_p)
  v__cl = ly.layerpotS(s=so).dot(phi)
  v__cl = submean(v__cl, so.w, 'Composite Linear')
  e = np.concatenate(( e, [test_N0_out(f_ex - v__cl, 'Composite Linear')] ))
  # --------------------------------------------------------------
  return e, f_ex, v__dp, v__d, v__clp, v__cl

def test_N0_total(n_ex):
  # f_ex, s_ex = exact(n_ex)
  rng = range(10 + (n_ex % 2), n_ex, 10)
  err = [np.zeros(4)]
  for k, n in enumerate(rng):
    print('')
    print('----------------------------------------------------------------------------------')
    print('------------------------------------ it %2d ---------------------------------------' %k)
    s = gm.sg_one_kite(n)
    result = test_N0(s)
    err = np.concatenate(( err, [np.array(result[0])] ))
  print('----------------------------------------------------------------------------------')
  if False:
    plt.figure()
    for k in range(4):
      plt.plot(err[:, k], '+-')
      plt.show(block=False)
      temp = input('Press')
  temp = input('Press')
  return rng, err


if __name__ == '__main__':
  n_ex = int(sys.argv[1])
  test_N0_total(n_ex)
