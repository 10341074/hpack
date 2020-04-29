from src import *
import dill


dill.load_session('dill3.pkl')
it_alpha = 15
p.ipb_opt_append(it_alpha=it_alpha)
p.alpha_fixed_ratio(0)
p.plot()
plt.savefig('prova32_0.svg')
p.alpha_fixed_ratio(-1)
p.plot()
plt.savefig('prova32_end.svg')

dill.dump_session('dill3-2.pkl')
# dill.load_session('dill.pkl')

end = input('Press enter')
# dill.load_session('dill.pkl')
