import plot as p
import numpy as np
import segment as s

r = 0.2
seg = s.Segment(n=100, f=s.circle, inargs=(0, r))
b = s.Boundary([seg])

l = s.Layer(b, s.ly.layerpotS)

(x, y, pp) = p.meshgrid((-2, 2, 100))
pp = s.Pointset(pp)

vv = 1.0 / 2 / s.pi / (r) * s.eval_layer(l, pp)
p.plot(x, y, vv)


pp = s.Pointset(np.arange(0.2, 0.2, 0.05))
vv = 1.0 / 2 / s.pi / (r) * s.eval_layer(l, pp)
