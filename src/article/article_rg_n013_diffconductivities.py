from src import *
import dill
# -------------------------
savefig = True
nd = 200
# filenames = ['two_ell_twoc']
filenames = ['two_kite_twoc']
# -----------------------------------------------------
hlist = [(2.0, 2.0), (2.0, 10.0), (2.0, 50.0), (2.0, 100.0), (100.0, 100.0)]
dpb.two_conductivities_flag = True
# -----------------------------------------------------
# dirtheta = ['0', '-p025pi', '-m025pi']
dirtheta = ['0']
# values_theta = [0, np.pi/4, -np.pi/4]
values_theta = [0]
# -----------------------------------------------------
alpha = 1e-12
h = 2.0

rb = 3.5
method='rg'
# -------------------------
for k_theta, theta in enumerate(values_theta):
  for filename in filenames:
    for hitem in hlist:
      p = m.EIT(theta=theta, alpha=alpha, h=h)
      p.domain(nsb=200, nso=200, rb=rb)
      p.ld = gm.template_lay_ld(filename, nd, h1=hitem[0], h2=hitem[1])
      # p.flag_points_inside(gm.template_lay_ld_flag(filename))

      if method == 'ls':
        p.solver()
      if method == 'rg':
        p.rg_solver()

      p.ipb()
      p.plot(pltype='c', levels=10, colorbar=False, linewidths=0.8, cmap='summer')
      plt.axis([-3, 3, -3, 3])
      plt.axis('square')
      plt.show(block=False)
      # ---- saving ----
      if savefig:
        plt.savefig('runs/%s_%s_theta%s_alpha%.0e_h1%d_h2%d.eps' %(method, filename, dirtheta[k_theta], alpha, hitem[0], hitem[1]), bbox_inches='tight')
      plt.close()
# ----------------
end = input('press enter')
