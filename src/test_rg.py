from src import *
p = m.EIT()
p.domain('one_kite')
p.meshgrid()
p.theta = 0
p.rg_solver()
p.rg_ipb()
p.plot()
inp = input('Press')


p.solver()
p.ipb()
p.plot()
inp = input('Press')
