from src import *
import tools_plots
savefig = False
plotstats = False
def funcnoise_lsm(p, noiselevelK, prename = (), name = ()):
  #########################################################################################
  p.noiselevelK = noiselevelK
  p.solver()
  p.ipb()  
  ###########################
  fig = plt.figure()
  plt.contour(p.x, p.y, p.z.reshape(p.y.size, p.x.size), 30)
  # p.ld.plot(lw = 1.2, ms = 1)
  p.ld.plot(lw = 2.8, ms = 1, color='m')
  p.so.plot(lw = 1.2, ms = 1)
  plt.axis('equal')
  plt.axis('square')
  plt.show(block=False)
  if savefig:
    plt.savefig('runs/fig-thesis/%s_lsm_%s_noisel%02.0e_noiselK%02.0e.eps' %(prename, name, p.noiselevel, p.noiselevelK), bbox_inches='tight')
  return
######################################################################################
# FACTORIZATION
######################################################################################
def funcnoise_fm(p, noiselevelK, m0 = 30, prename = (), name = ()):
  p.m0 = m0
  p.noiselevelK = noiselevelK
  p.fact_solver()
  p.LLfact = p.LL0
  p.fact_addnoiseK()
  p.fact_ieig()
  ######################################
  z = np.array(p.z)
  for k in range(p.x.size):
    if p.p.flag_inside_s[k] == 0:
      z[k] = 0
  z = z.reshape((p.x.size, p.y.size))
  fig = plt.figure()
  fig = plt.contourf(p.x, p.y, z, 1)
  # p.ld.plot(lw = 1.2, ms = 1)
  p.ld.plot(lw = 2.8, ms = 1, color='m')
  plt.axis('equal')
  plt.axis('square')
  plt.show(block=False)
  if savefig:
    plt.savefig('runs/fig-thesis/%s_fm_%s_m0%s_noisel%02.0e_noiselK%02.0e.eps' %(prename, name, p.m0, p.noiselevel, p.noiselevelK), bbox_inches='tight')
  #############################################
  # eigplot(p.fact_wsorted, p.m0, p.fact_linreg, prename = prename, name = name, nlK= p.noiselevelK)
  return
def funcnoise_fm_eig(p, noiselevelK, m0 = 30, prename = (), name = ()):
  fig = plt.figure()
  p.noiselevelK = 0
  ##########################################
  p.fact_solver()
  p.LLfact = p.LL0
  p.fact_addnoiseK()
  p.fact_ieig()
  eigplot(p.fact_wsorted, 0, p.fact_linreg, prename = prename, name = name, nlK= p.noiselevelK, marker = 'x', color = 'k', vec = [29])
  ###########################################
  p.m0 = m0
  p.noiselevelK = noiselevelK
  ##########################################
  p.fact_solver()
  p.LLfact = p.LL0
  p.fact_addnoiseK()
  p.fact_ieig()
  eigplot(p.fact_wsorted, p.m0, p.fact_linreg, prename = prename, name = name, nlK= p.noiselevelK, color = 'b', vec=[9])
  ###########################################
  return
 
def eigplot(wsorted, m0, linreg, prename = (), name = (), nlK = (), marker = 'o', color = 'k', vec = []):
  # plot
  x = np.arange(m0)
  if m0 > 0:
    plt.plot(x, linreg.intercept + linreg.slope*x, 'b-', ms=0.2, lw = 2)
  for ind in vec:
    temp_y = np.array([1, 1]) * np.log(wsorted[0][int(ind)])
    plt.plot([0, len(wsorted[0])], temp_y, 'k:', ms=0.2, color=color)
  # temp_y = np.array([1, 1]) * np.log(wsorted[0][int(29)])
  # plt.plot([0, len(wsorted[0])], temp_y, 'k:', ms=0.2, color=color)
  # temp_y = np.array([1, 1]) * np.log(wsorted[0][int(39)])
  # plt.plot([0, len(wsorted[0])], temp_y, 'k:', ms=0.2)
  plt.plot(range(len(wsorted[0])), np.log(wsorted[0]), marker, ms=5, markeredgewidth=1, markeredgecolor=color,
           markerfacecolor='None')
  plt.show(block=False)  
  if savefig:
    plt.savefig('runs/fig-thesis/%s_fm_%s_m0%s_%s_noiselK%02.0e.eps' %(prename, name, m0, "eig", nlK), bbox_inches='tight')
  return

if __name__ == '__main__':
  prename = 'noise'
  name = 'ellipse'
  p = m.EIT()
  # p.domain('three_ellipse', nsd = 100)
  p.domain()
  p.meshgrid((-3, 3, 80))
  p.noiselevel = 0
  # p.delta = p.noiselevel
  p.alpha = 1e-10
  funcnoise_lsm(p, prename = prename, name = name, noiselevelK = 1e-2)
  funcnoise_lsm(p, prename = prename, name = name, noiselevelK = 1e-3)
  funcnoise_fm(p, prename = prename, name = name, noiselevelK = 1e-3, m0 = 10)
  funcnoise_fm(p, prename = prename, name = name, noiselevelK = 1e-3, m0 = 12)
  funcnoise_fm(p, prename = prename, name = name, noiselevelK = 1e-3, m0 = 14)
  funcnoise_fm_eig(p, noiselevelK = 1e-3, m0 = 12, prename = prename, name = name)  
  ret = input('Press')
