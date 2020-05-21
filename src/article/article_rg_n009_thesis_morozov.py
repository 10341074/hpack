from src import *
import dill
# -------------------------
dump = True
filenames = ['kite_025pi_small', 'ell_thin', 'one_tri', 'two_ell', 'three_ell', 'one_tri_center']
dirtheta = ['theta0', 'theta-p025pi', 'theta-m025pi']
nd = 200
it_alpha = 1
alpha = 1e-14
h = 2.0
rb=3.8
# -------------------------
for k_theta, theta in enumerate([0, np.pi/4, -np.pi/4]):
  for filename in filenames:
    p = m.EIT(theta=theta, alpha=alpha, h=h)
    p.domain(nsb=200, nso=200, rb=rb)
    p.ld = gm.template_lay_ld(filename, nd)
    p.flag_points_inside(gm.template_lay_ld_flag(filename))

    p.rg_solver()

    p.ipb_opt(it_alpha=it_alpha)
    # ---- saving ----
    if dump:
      dill.dump_session('runs/' + dirtheta[k_theta] + '/' + '200502_morozov_' + filename + '.pkl')
# ----------------
end = input('press enter')
