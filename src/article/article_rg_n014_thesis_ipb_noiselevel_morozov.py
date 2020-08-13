from src import *
import dill
# -------------------------
dump = False
savefig = True
nd = 200
filenames = ['kite_025pi_small', 'ell_thin', 'one_tri', 'two_ell', 'three_ell', 'one_tri_center', 'one_square', 'one_cross', 'one_tri_square', 'one_mix']
# -----------------------------------------------------
# dirtheta = ['0', '-p025pi', '-m025pi']
dirtheta = ['-p025pi']
# values_theta = [0, np.pi/4, -np.pi/4]
values_theta = [np.pi/4]
# -----------------------------------------------------
alpha = [1e-12, 4.6e-13, 2.1e-13, 1e-13, 4.6e-14, 2.1e-14, 1e-14, 4.6e-15, 2.1e-15, 1e-15, 4.6e-16, 2.1e-16, 1e-16, 4.6e-17, 2.1e-17, 1e-17]
delta = 5e-8
delta = 5e-7
alpha = [1e-8, 4.6e-9, 2.1e-9, 1e-9, 4.6e-10, 2.1e-10, 1e-10, 4.6e-11, 2.1e-11, 1e-11, 4.6e-12, 2.1e-12, 1e-12, 4.6e-13, 2.1e-13, 1e-13]
delta = 1e-6
h = 2.0

rb = 3.5
method='rg'

noiselevel=0.0
noiselevelU=0.0
# ipb.sasaving_disc = 0
# -------------------------
for k_theta, theta in enumerate(values_theta):
  for filename in filenames:
    p = m.EIT(theta=theta, alpha=1e-12, h=h, noiselevel=noiselevel, noiselevelU=noiselevelU, delta=delta)
    p.domain(nsb=200, nso=200, rb=rb)
    p.ld = gm.template_lay_ld(filename, nd)
    p.flag_points_inside(p.so)
    # p.flag_points_inside(gm.template_lay_ld_flag(filename))

    if method == 'ls':
      p.solver()
    if method == 'rg':
      p.rg_solver()

    discrepancy = np.zeros((len(alpha), len(p.testpoints.x)))
    ratio       = np.zeros((len(alpha), len(p.testpoints.x)))
    for k_alpha, alpha_new in enumerate(alpha):
      p.solver_update_alpha(alpha_new)
      p.ipb()
      res = (p.isolver.BY.dot(p.isolver.A.dot(p.isolver.save_zeta.T) - p.isolver.save_rhs.T)).T
      discrepancy[k_alpha, :] = ipb.func_disc_value_pow1_n(res, p.isolver.save_sol, p.isolver.w, p.isolver.delta)
      ratio[k_alpha, :] = p.isolver.save_ratio[:, 0]
    z = ratio[np.argmin(abs(discrepancy), axis=0), range(ratio.shape[1])]
    # --------------------------------------------
    filestring = 'runs/%s_%s_theta%s_alpha%.0e_h%d' %(method, filename, dirtheta[k_theta], alpha_new, h)
    if savefig:
      p.plot(z, pltype='c', levels=10, colorbar=False, linewidths=0.8, cmap='summer')
      plt.axis([-3, 3, -3, 3])
      plt.axis('square')
      plt.show(block=False)
      # ---- saving ----
      plt.savefig(filestring + '.eps', bbox_inches='tight')
      plt.close()
      # ----------------
      p.plot(np.argmin(abs(discrepancy), axis=0), vmin=0, vmax=len(alpha)-1)
      plt.savefig(filestring + '_index' + '.eps', bbox_inches='tight')
      plt.close()
    if dump:
      dill.dump_session(filestring + '.pkl')
# ----------------
end = input('press enter')
