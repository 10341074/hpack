from src import *
import dill
# -------------------------
savefig = True
nd = 200
filenames = ['kite_025pi_small', 'ell_thin', 'one_tri', 'two_ell', 'three_ell', 'one_tri_center', 'one_square', 'one_cross', 'one_tri_square', 'one_mix']
# filenames = ['two_kite', 'one_kite', 'kite_025pi_small']
# -----------------------------------------------------
# dirtheta = ['0', '-p025pi', '-m025pi']
dirtheta = ['0']
# values_theta = [0, np.pi/4, -np.pi/4]
values_theta = [0]
# -----------------------------------------------------
alpha = [1e-8, 1e-12, 1e-14]
h = 2.0

rb = 3.5
method='rg'

noiselevel=0.0
noiselevelU=0.03
# -------------------------
for k_theta, theta in enumerate(values_theta):
  for filename in filenames:
    p = m.EIT(theta=theta, alpha=1e-12, h=h, noiselevel=noiselevel, noiselevelU=noiselevelU)
    p.domain(nsb=200, nso=200, rb=rb)
    p.ld = gm.template_lay_ld(filename, nd)
    # p.flag_points_inside(gm.template_lay_ld_flag(filename))

    if method == 'ls':
      p.solver()
    if method == 'rg':
      p.rg_solver()

    for alpha_new in alpha:
      p.solver_update_alpha(alpha_new)
      p.ipb()
      p.plot(pltype='c', levels=10, colorbar=False, linewidths=0.8, cmap='summer')
      plt.axis([-3, 3, -3, 3])
      plt.axis('square')
      plt.show(block=False)
      # ---- saving ----
      if savefig:
        plt.savefig('runs/%s_%s_theta%s_alpha%.0e_h%d.eps' %(method, filename, dirtheta[k_theta], alpha_new, h), bbox_inches='tight')
      plt.close()
# ----------------
end = input('press enter')
