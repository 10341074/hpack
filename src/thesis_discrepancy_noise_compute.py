from src import *
import dill

dump = True

p = m.EIT()
p.domain()
# p.ld = gm.two_ellipse(100)
p.ld = gm.three_ellipse(100)
# nd = 100; p.ld = sg.Layer([sg.Boundary([sg.Segment(nd, f_inargs = (sh.ellipse, (0, 0.1, 1.2)), quad='ps', aff=(-1, np.exp(1j * np.pi * 2/3))) ])])



p.meshgrid((-3, 3, 80))
p.noiselevel = 0
p.noiselevelK = 1e-3 ###################### NOISE K
p.alpha = 1e-10
p.delta = 1e-6
p.solver()

it_alpha = 20
p.ipb_opt(it_alpha=it_alpha)

if dump:
  # dill.dump_session('runs/dill/dill_two_ellipse_delta1e-6.pkl')
  # dill.dump_session('runs/dill/dill_three_ellipse_delta1e-6.pkl')
  # dill.dump_session('runs/dill/dill_sharp_ellipse_delta1e-6.pkl')
  dill.dump_session('runs/dill/dill_three_ellipse_noiseK1e-3_delta1e-6_bis.pkl')
dill.dump_session('prova.pkl')

end = input('Press enter')
