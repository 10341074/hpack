from src import *
import dill
# -------------------------
dump = False
filename = 'prova'
nd = 100
it_alpha = 2
# -------------------------

p = m.EIT()
p.domain()
# p.ld = gm.two_ellipse(nd)
# p.ld = gm.three_ellipse(nd)
p.ld = sg.Layer([sg.Boundary([sg.Segment(nd, f_inargs = (sh.ellipse, (0, 0.1, 1.2)), quad='ps', aff=(-1, np.exp(1j * np.pi * 2/3))) ])])

p.meshgrid((-3, 3, 80))
p.noiselevel = 0
# p.noiselevelK = 1e-3 # ''' noise K '''
p.alpha = 1e-10
p.delta = 1e-6
p.rg_solver()

p.ipb_opt(it_alpha=it_alpha)

# ---- saving ----
if dump:
  dill.dump_session('runs/dill/' + filename + '.pkl')
dill.dump_session(filename + '.pkl')
# ----------------
end = input('press enter')
