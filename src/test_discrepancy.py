from src import *
import dill

p = m.EIT()
p.domain()
p.ld =gm.three_ellipse(100)
 
p.meshgrid((-2, 2, 40))
p.solver()
# alpha = np.array([1e-16 + 1e-16 * k for k in range(10)])
# alpha = np.concatenate((alpha, [1e-15*10**k for k in range(8)]))
#p.test_alpha(alpha=alpha)

#
it_alpha = 15
p.ipb_opt(it_alpha=it_alpha)
p.alpha_fixed_ratio(0)
p.plot()
# plt.savefig('prova3_0.svg')
p.alpha_fixed_ratio(-1)
p.plot()
# plt.savefig('prova3_end.svg')

# dill.dump_session('dill3.pkl')
# dill.load_session('dill.pkl')

end = input('Press enter')
