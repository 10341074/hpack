from src import *
import dill

dump = False

p = m.EIT()
p.domain()
p.ld = gm.three_ellipse(100)

p.meshgrid((-3, 3, 80))
# p.noiselevel = 0
# p.alpha = 1e-10
# p.delta = 1e-6
p.solver()

it_alpha = 3
p.ipb_opt(it_alpha=it_alpha)

if dump:
  dill.dump_session('runs/dill/dill_three_ellipse_delta1e-6.pkl')
dill.dump_session('prova.pkl')

end = input('Press enter')
