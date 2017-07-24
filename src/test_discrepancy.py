from src import *
import dill

p = m.EIT()
p.domain()
p.ld = gm.two_ellipse(40)

p.meshgrid((-2, 2, 20))
p.noiselevel = 0
p.alpha = 1e-10
p.delta = 1e-6
p.solver()

it_alpha = 4
p.ipb_opt(it_alpha=it_alpha)
p.alpha_fixed_ratio(-1)

plt.contour(p.x, p.y, p.z.reshape(p.y.size, p.x.size), 30)
plt.axis('equal')
plt.axis('square')
plt.show()
dill.dump_session('runs/dill/dill_two_ellipse_delta1e-6.pkl')
end = input('Press enter')
