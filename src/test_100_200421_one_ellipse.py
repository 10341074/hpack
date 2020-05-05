from src import *
import tools_plots

p = m.EIT()
p.domain('one_ellipse', nsb=100, nso=80, nsd=40)
# p.ld = sg.Layer([sg.Boundary([sg.Segment(40, f_inargs = (sh.ellipse, (0, 2, 1)), quad='ps') ])])
# p.ld = gm.sg_one_triangle(40)
p.meshgrid((-2, 2, 20))
# p.so.BX = linf.base_mean(p.so.BX, p.so.w)
print('--------------------------------------------------------')
print('--------------------- solver ---------------------------')
p.solver()
print('--------------------------------------------------------')
print('--------------------- ipb ------------------------------')
p.ipb()
print('--------------------------------------------------------')
print('--------------------- plot -----------------------------')
p.plot()
# -----------------------------------------------------------------------------------------------


# tools_plots.plot_contourf(p, ())

# alpha = np.array([1e-16 + 1e-16 * k for k in range(10)])
# alpha = np.concatenate((alpha, [1e-15*10**k for k in range(8)]))
#p.test_alpha(alpha=alpha)



# L0 = computeL0(p.so, p.so.BX, p.pp)
# z = L0.dot(p.isolver.save_sol[0])
# p.plot(z)

###########################
# p.solver()
# it_alpha = 8
# p.ipb_opt(it_alpha=it_alpha)
# p.alpha_fixed_ratio(0)
# p.plot()
# p.alpha_fixed_ratio(-1)
# p.plot()
#############################
# dill.dump_session('dill.pkl')
# dill.load_session('dill.pkl')

# CORRECT
###########################
# p.fact_solver()
# p.fact_ieig()
# v = tools_plots.plot_contourf_1(p, ())
# p.plot()
# plt.plot(v[:,0], v[:,1],'r-')
# o = tools_plots.stats(v)
# plt.plot([o[0].real, v[o[3],0]], [o[0].imag, v[o[3],1]], 'b-')
# plt.show(block=False)

end = input('Press enter')
