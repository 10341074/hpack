from src import *
''' This script checks computeU in source code is correct, and compares two methods, with y1, y2 evaluated on so '''
p = m.EIT(h=200)
p.domain(nsb=4, nso=150, nsd=150, rb=3.5)
p.ld = gm.template_lay_ld('one_tri', 150)
# ---------------------------------------------------
''' method 1: use 'computeU' from module ''' 
allrhs, allpsi = ipb.computeallpsi(p.ld, p.sb, p.c)
U, U_nu = ipb.computeU(allpsi, p.ld, p.so, p.sb)
# ---------------------------------------------------
''' method 2: manual compute U with allpsi'''
def func(targets = p.testpoints):
    return (ly.layerpotS(s=p.ld, t=targets).dot(allpsi) + ly.layerpotS(s=p.sb, t=targets, noweights=True))[:, 0]
z = func()
''' plot u on 2d meshgrid '''
p.plot(z, pltype='cf', levels=(np.log(np.linspace(1, 2, 40)) / np.log(2)) * (max(z) - min(z)) + min(z))
# ---------------------------------------------------
plt.figure()
y2 = np.diag(p.so.w).dot(func(p.so))
y1 = U[0, :]
plt.plot(y1, '1-')
plt.plot(y2, '2-')
plt.show(block=False)