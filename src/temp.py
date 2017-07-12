from src import *
import test_plots
p = m.EIT()
p.domain()
p.meshgrid()
p.plot_domain()
p.noiselevel = 1e-2
p.delta = p.noiselevel
p.solver()
p.ipb()
p.plot()

prename = 'noise'
name = 'ellipse'
k = 0

fig = plt.figure()
plt.contour(p.x, p.y, p.z.reshape(p.y.size, p.x.size), 30)
p.ld.plot(lw = 0.8, ms = 1)
p.so.plot(lw = 0.8, ms = 1)
plt.axis('equal')
plt.axis('square')
# plt.show(block=False)
plt.savefig('runs/fig-thesis/%s_lsm_%s%s_noiselevel%02.0e.eps' %(prename, name, k, p.noiselevel), bbox_inches='tight')
#%02.03e
plt.figure()
p.fact_solver()
p.fact_ieig()

p.plot()

# v = test_plots.plot_contourf_1(p, ())
# #plt.plot(v[:,0], v[:,1],'r-')
# o = test_plots.stats(v)
# #plt.plot([o[0].real, v[o[3],0]], [o[0].imag, v[o[3],1]], 'b-')
# print(o)
# p.ld.plot(lw = 0.8, ms = 1)
# plt.axis('equal')
# plt.axis('square')

# plt.savefig('runs/fig-thesis/%s_fm_%s%s_noiselevel%s.eps' %(prename, name, k, p.noiselevel), bbox_inches='tight')


i = input('Press')
