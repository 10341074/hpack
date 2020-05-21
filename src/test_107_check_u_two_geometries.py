from src import *
''' Guess: u, u_nu depends on the inclusion, but the mixed integral at the rhs in rg: (u * f_nu - u_nu * f), maybe not '''
p_so = m.EIT()
p_so.domain(nso=400)
so = p_so.so
# ---------------------------------------
z0 = 0
theta = 0.2
# ----------
F    = (np.cos(theta) * ly.phi_x(z0, so.x)          + np.sin(theta) * ly.phi_y(z0, so.x))
F_nu = (np.cos(theta) * ly.phi_x_n(z0, so.x, so.nx) + np.sin(theta) * ly.phi_y_n(z0, so.x, so.nx))
# ---------------------------------------
def compU(geom, h=2):
    p = m.EIT(h=h)
    p.domain(nsb=1, nso=400, nsd=150, rb=3.5, quad='p')
    p.ld = gm.template_lay_ld(geom, 150)

    allrhs, allpsi = ipb.computeallpsi(p.ld, p.sb, p.c)
    U, U_nu = ipb.computeU(allpsi, p.ld, p.so, p.sb)
    return U, U_nu


u1, unu1 = compU('one_ell', h=100)
u2, unu2 = compU('one_ell')
u3, unu3 = compU('one_tri')
plt.figure()
plt.plot(u1[0], '1-')
plt.plot(u2[0], '1-')
plt.plot(u3[0], '2-')
plt.show(block=False)

print(u1.dot(F_nu) - unu1.dot(F))
print(u2.dot(F_nu) - unu2.dot(F))
print(u3.dot(F_nu) - unu3.dot(F))
