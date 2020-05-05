from src import *
import tools_plots

savefig = False
# ======================================================================================================================
# --- initialization ---
def init_geomanal_geomnum_std():
  p = m.EIT()
  p.domain()
  p.meshgrid((-2, 2, 60))
  p.sb = sg.Segment(80, f_inargs = (sh.circle, (0, 6)), quad='ps')
  # ----------------------------
  # --- settings ---
  p.alpha = 1e-15
  p.theta = 0
  return p
# ======================================================================================================================
def compare_ls_vs_rg(p, rg = True, ls = True):
  '''
    with this function we want to compare: linearsampling vs rg
  '''
  # ----------------------------
  if rg:
    p.rg_solver()
    p.ipb()
    p.plot()
    plt.title('reciprocal gap'); plt.show(block=False)
  # ==========================
  if ls:
    p.solver()
    p.ipb()
    p.plot()
    plt.title('linear sampling'); plt.show(block=False)
  return p
# ====================================================================================================================
if __name__ == "__main__":
  if True:
    p = init_geomanal_geomnum_std()
    compare_ls_vs_rg(p)
  inp = input('press')

