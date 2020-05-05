from src import *
import dill
# -------------------------
dump = False
filename = 'prova_alphafixed'
nd = 100
it_alpha = 2
# -------------------------

p = m.EIT()
p.domain(nsb=200, nso=200)
p.ld = gm.template_lay_ld('ell_thin', 200)

p.rg_solver()

p.ipb_alpha_fixedarray()

# ---- saving ----
if dump:
  dill.dump_session('runs/dill/' + filename + '.pkl')
dill.dump_session(filename + '.pkl')
# ----------------
end = input('press enter')
