from src import *
import test_plots

savefig = False
plotstats = False
my_fact_L_trigbasis = False



ellipses = [
    sg.Layer([sg.Boundary([sg.Segment(80, f_inargs = (sh.ellipse, (0, 2, 1)), quad='ps') ])]), 
    sg.Layer([sg.Boundary([sg.Segment(80, f_inargs = (sh.ellipse, (0, 1.8, 1.2)), quad='ps', aff=(1, np.exp(1j * np.pi/4))) ])]),
    sg.Layer([sg.Boundary([sg.Segment(80, f_inargs = (sh.ellipse, (0, 0.1, 1.2)), quad='ps', aff=(-1, np.exp(1j * np.pi * 2/3))) ])])
]

def first_compare(ld = (), prename = (), name = () ):
  if ld == ():
    n = len(ellipses)
    prename, name = 'one_ellipse', 'ellipse'
  else:
    n = 1
  for k in range(0, n):
    p = m.EIT()
    p.domain('one_ellipse', nso = 150)
    # p.ld = sg.Layer([sg.Boundary([sg.Segment(40, f_inargs = (sh.ellipse, (0, 2, 1)), quad='ps') ])])
    # p.ld = sg.Layer([sg.Boundary([sg.Segment(80, f_inargs = (sh.ellipse, (0, 1.5, 1)), quad='ps', aff=(1, 1.2*np.exp(1j * np.pi/4))) ])])
    if ld == ():
      p.ld = ellipses[k]
    else:
      p.ld = ld
    p.meshgrid((-3, 3, 120))
    ############################################
    p.noiselevel = 0
    #############################################
    # FACTORIZATION METHOD
    ##############################################
    # if my_fact_L_trigbasis:
    #   setups.fact_L_trigbasis = True
    p.m0 = 30
    p.fact_solver()
    p.LLfact = p.LL0
    print("Trig is", setups.fact_L_trigbasis)
    p.fact_ieig()
    #############################################
    # plot
    v = test_plots.plot_contourf_1(p, ())
    o = test_plots.stats(v)
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
    ############################################################################
    ############################################################################
    # # LLdiff
    # p.m0 = 30
    # p.fact_solver()
    # p.LLfact = p.LLdiff
    # print("Trig is", setups.fact_L_trigbasis)
    # p.fact_ieig()
    # #############################################
    # # plot
    # v = test_plots.plot_contourf_1(p, ())
    # o = test_plots.stats(v)
    # if plotstats:
    #   plt.plot(v[:,0], v[:,1],'b-')
    #   # plt.plot([o[0].real, v[o[3],0]], [o[0].imag, v[o[3],1]], 'b-')
    # print(o[0].real, o[0].imag, o[1], o[2])
    # p.ld.plot(lw = 1.2, ms = 1)
    # plt.axis('equal')
    # plt.axis('square')
    # plt.show(block=False)
    # if savefig:
    #   plt.savefig('runs/fig-thesis/%s_fm_%s%s_LLdiff.eps' %(prename, name, k), bbox_inches='tight')
    # #############################################
    # eigplot(p.fact_wsorted, p.m0, p.fact_linreg)
    ################################################################################
    p.m0 = 40
    p.fact_solver()
    p.LLfact = p.LL0
    print("Trig is", setups.fact_L_trigbasis)
    p.fact_ieig()
    #############################################
    # plot
    v = test_plots.plot_contourf_1(p, ())
    o = test_plots.stats(v)
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
    ###########################################
    eigplot(p.fact_wsorted, p.m0, p.fact_linreg, prename = prename, name = name)
  return

def eigplot(wsorted, m0, linreg, prename = (), name = ()):
  # plot
  fig = plt.figure()
  x = np.arange(m0)
  plt.plot(x, linreg.intercept + linreg.slope*x, 'b-', ms=0.2)
  temp_y = np.array([1, 1]) * np.log(wsorted[0][int(29)])
  plt.plot([0, len(wsorted[0])], temp_y, 'k:', ms=0.2)
  temp_y = np.array([1, 1]) * np.log(wsorted[0][int(39)])
  plt.plot([0, len(wsorted[0])], temp_y, 'k:', ms=0.2)
  plt.plot(range(len(wsorted[0])), np.log(wsorted[0]), 'ko', ms=2, markeredgewidth=1, markeredgecolor='k',
           markerfacecolor='None')
  plt.show(block=False)  
  if savefig:
    plt.savefig('runs/fig-thesis/%s_fm_%s_m0%s_%s.eps' %(prename, name, m0, "eig"), bbox_inches='tight')
  return


if __name__ == "__main__":
  first_compare(gm.two_ellipse(100), 'first_compare', 'two_ellipse')
  first_compare(gm.three_ellipse(100), 'first_compare', 'three_ellipse')
  
  end = input('Press enter')
