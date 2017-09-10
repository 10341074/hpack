from src import *
import test_plots
savefig = False
plotstats = False
def funcnoise():
  p = m.EIT()
  p.domain()
  p.meshgrid((-3, 3, 80))
  p.noiselevel = 1e-2
  p.delta = p.noiselevel
  p.alpha = 1e-10
  p.solver()
  p.ipb()
  prename = 'noise'
  name = 'ellipse'
  
  ###########################
  fig = plt.figure()
  plt.contour(p.x, p.y, p.z.reshape(p.y.size, p.x.size), 30)
  p.ld.plot(lw = 1.2, ms = 1)
  p.so.plot(lw = 1.2, ms = 1)
  plt.axis('equal')
  plt.axis('square')
  plt.show(block=False)
  if savefig:
    plt.savefig('runs/fig-thesis/%s_lsm_%s%s_noiselevel%02.0e.eps' %(prename, name, k, p.noiselevel), bbox_inches='tight')
  #
  ######################################################################################
  # FACTORIZATION
  ######################################################################################
  p.m0 = 30
  p.fact_solver()
  p.LLfact = p.LL0
  p.fact_ieig()
  ######################################
  z = np.array(p.z)
  for k in range(p.x.size):
    if p.p.flag_inside_s[k] == 0:
      z[k] = 0
  z = z.reshape((p.x.size, p.y.size))
  fig = plt.figure()
  fig = plt.contourf(p.x, p.y, z, 1)

  p.ld.plot(lw = 1.2, ms = 1)
  plt.axis('equal')
  plt.axis('square')
  plt.show(block=False)
  if savefig:
    plt.savefig('runs/fig-thesis/%s_fm_%s%s.eps' %(prename, name, k), bbox_inches='tight')
    #############################################
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

if __name__ == '__main__':
  funcnoise()
  ret = input('Press')
