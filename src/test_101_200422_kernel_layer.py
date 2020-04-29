from src import *
import layerpot_test as ly_t

s = sg.Segment(400, f_inargs = (sh.circle, (0, 3)), quad='ps')
# s = sg.Segment(400, f_inargs = (sh.ellipse, (0, 3, 4)), quad='ps')
ly_t.plot_kernel_SegmentPoint(s, layertype='Kp')
plt.show(block=True)