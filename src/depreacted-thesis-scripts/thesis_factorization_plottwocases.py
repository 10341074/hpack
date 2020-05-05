from src import *
import tools_plots

savefig = False
plotstats = False
my_fact_L_trigbasis = False

def first_compare(ld = (), prename = (), name = () ):
  p = m.EIT()
  p.domain('one_ellipse', nso = 150)
  p.ld = ld
  p.meshgrid((-3, 3, 120))
  ############################################
  p.noiselevel = 0
  #############################################
  # FACTORIZATION METHOD
  ##############################################
  # if my_fact_L_trigbasis:
  #   setups.fact_L_trigbasis = True
  p.m0 = 40
  p.fact_solver()
  p.LLfact = p.LL0
  print("Trig is", setups.fact_L_trigbasis)
  p.fact_ieig()
  #############################################
  # plot
  v = tools_plots.plot_contourf_1(p, ())
  o = tools_plots.stats(v)
  if plotstats:
    plt.plot(v[:,0], v[:,1],'b-')
    # plt.plot([o[0].real, v[o[3],0]], [o[0].imag, v[o[3],1]], 'b-')
  print(o[0].real, o[0].imag, o[1], o[2])
  p.ld.plot(lw = 1.2, ms = 1)
  plt.axis('equal')
  plt.axis('square')
  plt.show(block=False)
  if savefig:
    plt.savefig('runs/fig-thesis/%s_fm_%s%s_m0%s.eps' %(prename, name, k, p.m0), bbox_inches='tight')
  #############################################
  eigplot(p.fact_wsorted, p.m0, p.fact_linreg, prename = prename, name = name)
  ############################################
  w, v, wsorted, m0, linreg = ipb.eigselect(p.LLfact, m0=p.m0)
  test_z, test_ws, test_linreg = ipb.testieig(w, v, wsorted, m0, linreg, p.isolver, sg.Pointset([-1.5+2j]), p.LL0, p.so, p.theta)
  print('Tested', test_z)
  y = test_ws - test_ws[0]
  plt.plot(range(len(test_ws)), y, 'k^', ms=6, markeredgewidth=1, markeredgecolor='k',
           markerfacecolor='None')
  x = np.arange(m0)
  plt.plot(x, test_linreg.intercept + test_linreg.slope*x  - test_ws[0], 'k-', ms=0.2)
  plt.show(block=False)
  ############################################
  w, v, wsorted, m0, linreg = ipb.eigselect(p.LLfact, m0=p.m0)
  test_z, test_ws, test_linreg = ipb.testieig(w, v, wsorted, m0, linreg, p.isolver, sg.Pointset([-1.23+0.03j]), p.LL0, p.so, p.theta)
  print('Tested', test_z)
  y = test_ws - test_ws[0]
  plt.plot(range(len(test_ws)), y, 'kv', ms=6, markeredgewidth=1, markeredgecolor='k',
           markerfacecolor='None')
  x = np.arange(m0)
  plt.plot(x, test_linreg.intercept + test_linreg.slope*x  - test_ws[0], 'k-', ms=0.2)
  plt.show(block=False)
  if savefig:
    plt.savefig('runs/fig-thesis/%s_fm_%s_m0%s_%s.eps' %(prename, name, m0, "eig"), bbox_inches='tight')
  return

def eigplot(wsorted, m0, linreg, prename = (), name = ()):
  # plot
  fig = plt.figure()
  x = np.arange(m0)
  plt.plot(x, linreg.intercept + linreg.slope*x- np.log(wsorted[0][0]), 'k-', ms=0.2)
  y = np.log(wsorted[0]) - np.log(wsorted[0][0])
  plt.plot(range(len(wsorted[0][:m0])), y[:m0], 'ko', ms=2.5, markeredgewidth=1, markeredgecolor='k',
           markerfacecolor='None')

  plt.show(block=False)
  if savefig:
    plt.savefig('runs/fig-thesis/%s_fm_%s_m0%s_%s.eps' %(prename, name, m0, "eig"), bbox_inches='tight')
  return


if __name__ == "__main__":
  first_compare(gm.one_ellipse(100), 'first_compare_factcases', 'one_ellipse')
  end = input('Press enter')
