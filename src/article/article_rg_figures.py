from src import *
savefig = False

def test_ipb(p=(), alpha=1e-5, nsb=80, name=""):
  if p == ():
    p = m.EIT()
  else:
    print("Imported p")
  p.domain()
  p.sb = sg.Segment(nsb, f_inargs = (sh.circle, (0, 6)), quad='ps')
  p.alpha = alpha
  p.meshgrid((-3, 3, 70))
  p.flag_points_inside(sg.Segment(200, f_inargs = (sh.circle, (0, 2.7)), quad='ps'))
  p.theta = 0
  p.rg_solver()
  p.ipb()
  p.plot()
  if savefig:
    plt.savefig('runs/fig-article/%s_rg_alpha_%s_sbn_%s_son_%s_ldn_%s.eps' %(name, p.alpha, p.sb.n, p.so.n, p.ld.n), bbox_inches='tight')
  p.solver()
  p.ipb()
  p.plot()
  if savefig:
    plt.savefig('runs/fig-article/%s_ls_alpha_%s_sbn_%s_son_%s_ldn_%s.eps' %(name, p.alpha, p.sb.n, p.so.n, p.ld.n), bbox_inches='tight')
  return p

def norms(K):
  c = numpy.linalg.cond(np.array(K, float))
  n2 = numpy.linalg.norm(np.array(K, float), ord=2)
  nI = numpy.linalg.norm(np.array(K, float), ord=np.inf)
  # print("Cond num= ", numpy.linalg.cond(np.array(K, float)), ' with for shape ', K.shape)
  # print("Norm 2=", numpy.linalg.norm(np.array(K, float), ord=2))
  # print("Norm inf=", numpy.linalg.norm(np.array(K, float), ord=np.inf))
  return c, n2, nI
def condnum_ipb(p=(), alpha=1e-5, nsb=100, name="twomethods", rg={}, ls={}):
  if p == ():
    p = m.EIT()
  else:
    print("Imported p")
  p.domain(nsb=nsb, nso=nsb)
  p.sb = sg.Segment(nsb, f_inargs = (sh.circle, (0, 6)), quad='ps')
  p.alpha = alpha
  # p.meshgrid((-2, 2, 60))
  p.theta = 0
  p.rg_solver()
  print("Cond num rg = ", numpy.linalg.cond(np.array(p.K, float)), ' with for shape ', p.K.shape)
  rg["c"], rg["n2"], rg["nI"] = norms(p.K)
  p.solver()
  print("Cond num ls = ", numpy.linalg.cond(np.array(p.K, float)), ' with for shape ', p.K.shape)
  norms(p.K)
  ls["c"], ls["n2"], ls["nI"] = norms(p.K)
  return p
def condnum_sens(p=(), name="condnum"):
  nvec = np.arange(40,190,6)
  # rg_cond, rg_norm2, rg_normI = np.empty(len(nvec)), np.empty(len(nvec)), np.empty(len(nvec))
  # ls_cond, ls_norm2, ls_normI = np.empty(len(nvec)), np.empty(len(nvec)), np.empty(len(nvec))
  rg = [{} for n in range(len(nvec))]
  ls = [{} for n in range(len(nvec))]
  for k, n in enumerate(nvec):
    condnum_ipb(nsb = nvec[k], rg = rg[k], ls = ls[k])
  plt.figure()
  plt.plot(nvec, [rg[k]["c"] for k in range(len(nvec))], '+-')
  plt.plot(nvec, [ls[k]["c"] for k in range(len(nvec))], '+-')
  plt.show(block=False)
  if savefig:
    plt.savefig('runs/fig-article/%s_condnum.eps' %(name), bbox_inches='tight')
                                                                             
  plt.figure()
  plt.plot(nvec, [rg[k]["n2"] for k in range(len(nvec))], '+-')
  plt.plot(nvec, [ls[k]["n2"] for k in range(len(nvec))], '+-')
  plt.title("Norm for p=2 of the discrete map")
  plt.show(block=False)
  if savefig:
    plt.savefig('runs/fig-article/%s_norm2.eps' %(name), bbox_inches='tight')

  plt.figure()
  plt.plot(nvec, [rg[k]["nI"] for k in range(len(nvec))], '+-')
  plt.plot(nvec, [ls[k]["nI"] for k in range(len(nvec))], '+-')
  plt.title("Norm for p=Inf of the discrete map")
  plt.show(block=False)
  if savefig:
    plt.savefig('runs/fig-article/%s_normI.eps' %(name), bbox_inches='tight')

  plt.figure()
  plt.plot(nvec, [ls[k]["c"] for k in range(len(nvec))], '+-')
  plt.show(block=False)
  if savefig:
    plt.savefig('runs/fig-article/%s_condnum_ls.eps' %(name), bbox_inches='tight')

if __name__ == "__main__":
  # test_ipb()
  # condnum_ipb()
  # condnum_sens()
  end = input('Press enter')
