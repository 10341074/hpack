import __test__ as t

from __load__ import *

s = sg.Segment(20, periodic=True, f=sg.circle, inargs=(-1, 0.4))
s2 = sg.Segment(20, periodic=True, f=sg.ellipse, inargs=(1, 0.4))
x, y, pp = plot.meshgrid((-3, 3, 60))
l = sg.Layer(b=[s], exp=ly.layerpotSD)
# l = sg.Layer(b=[s, s2], exp=ly.layerpotSD)
vv = sg.eval_layer(l, sg.Pointset(pp), exp=ly.layerpotS)
plot.plot(x, y, vv,t='srf')

#s.plot(p=True)
plt.show()
# A = ly.layerpotS(s=s)

#A = l.eval_self()


# >>> import time
# >>> print(time.time(), time.clock())
# 1359147652.31 0.021184
# >>> time.sleep(1)
# >>> print(time.time(), time.clock())
# 1359147653.31 0.02168

sd = sg.Segment(100, Z=sg.dZ, Zp=sg.dZp, Zpp=sg.dZpp, args=[], quad='gp')
sd.plot(p=True)
plt.show()
A = ly.layerpotD(s=sd)

print(A.dot(np.ones(len(A))))
