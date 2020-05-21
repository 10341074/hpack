from src import *
''' 
This module computes mean for F, and F_nu, increasing points in boundary
- the error decreases
- check if theory implies mean(F) = 0 in two dimension (the numeric integral seems to vanish)
'''
# -----------------------------------------------------------------------------------------------
theta = 0.2
z0 = 2.8j + 0.5
ro=3
# -----------------------------------------------------------------------------------------------
sf = []
sfnu = []

for k in np.arange(100, 300, 10):
    print(k)    
    so = sg.Segment(k, f_inargs = (sh.circle, (0, ro)), quad='ps')
    F    = np.cos(theta) * ly.phi_x(z0, so.x)          + np.sin(theta) * ly.phi_y(z0, so.x)
    F_nu = np.cos(theta) * ly.phi_x_n(z0, so.x, so.nx) + np.sin(theta) * ly.phi_y_n(z0, so.x, so.nx)
    sf   += [np.sum(F) / len(F)]
    sfnu += [np.sum(F_nu) / len(F_nu)]
plt.plot(np.abs(sf))
plt.yscale('log')
plt.show(block=False)