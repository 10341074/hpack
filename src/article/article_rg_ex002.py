from src import *
import tools_plots

savefig = False
# ----------------------------------------------------------------------------------------------------------------------
nb, no, nd = 300, 300, 200
quad = 'ps'
# -------------------------
ellipses = [
      # gm.lay_three_ellipse(nd)
      # sg.Layer([sg.Boundary([sg.poly([2, 1 +1.5j, -1 +1j], int(nd)) ])])
      sg.Layer([sg.Boundary([sg.Segment(nd, Z_args = (sh.kZ, sh.kZp, sh.kZpp, ()), quad='ps', aff=(0, np.exp(1j * np.pi / 4)))  ])])
      # sg.Layer([sg.Boundary([sg.Segment(nd, f_inargs = (sh.ellipse, (0, 0.1, 1.2)), quad='ps', aff=(0, np.exp(1j * np.pi * 2/3)))   ])])
      # sg.Layer([sg.Boundary([sg.Segment(nd, f_inargs = (sh.ellipse, (0, 2, 1)), quad='ps')                                          ])]) 
      # sg.Layer([sg.Boundary([sg.Segment(nd, f_inargs = (sh.ellipse, (0, 1.8, 1.2)), quad='ps', aff=(1, np.exp(1j * np.pi/4)))       ])])
      # sg.Layer([sg.Boundary([sg.Segment(nd, f_inargs = (sh.ellipse, (0, 0.1, 1.2)), quad='ps', aff=(0, np.exp(1j * np.pi * 2/3)))   ])])
      # sg.Layer([sg.Boundary([sg.Segment(nd, f_inargs = (sh.ellipse, (0, 0.1, 1.2)), quad='ps', aff=(-1, np.exp(1j * np.pi * 2/3)))  ])])
      # sg.Layer([sg.Boundary([sg.Segment(80, f_inargs = (sh.ellipse, (0, 1.5, 1)), quad='ps', aff=(1, 1.2*np.exp(1j * np.pi/4)))     ])])    
]
# ======================================================================================================================
# --- initialization ---
def init_geomanal_geomnum_2():
  p = m.EIT()
  p.domain('one_ellipse', nso = no, quad=quad)
  p.meshgrid((-3, 3, 80))
  # ----------------------------
  # --- settings ---
  p.alpha = 1e-12
  p.noiselevel = 0
  p.theta = np.pi / 4 
  # p.c = 100
  return p
# ======================================================================================================================
def run_rg(p, prename = (), name = () ):
  radiusb = 4
  n = len(ellipses)
  prename, name = 'one_ellipse', 'ellipse'
  for k in range(0, n):
    p.ld = ellipses[k]
    p.sb = sg.Segment(nb, f_inargs = (sh.circle, (0, radiusb)), quad=quad)
    # ----------------------------
    p.rg_solver()
    p.ipb()
    # ----------------------------
    # --- plot ----    
    fig = plt.figure()
    plt.contour(p.x, p.y, p.z.reshape(p.y.size, p.x.size), 30)
    p.ld.plot(lw = 1.2, ms = 1)
    p.so.plot(lw = 1.2, ms = 1)
    plt.axis('equal')
    plt.axis('square')
    plt.show(block=False)
    if savefig:
      plt.savefig('runs/fig-thesis/%s_rg_%s%s_alpha%s_no%snd%snb%s_radius%s.eps' %(prename, name, k, p.alpha, no, nd, p.sb.n, radiusb), bbox_inches='tight')
  return
# ====================================================================================================================
if __name__ == "__main__":
  if True:
    p = init_geomanal_geomnum_2()
    run_rg(p)
  inp = input('press')

