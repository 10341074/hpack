import numpy as np
import matplotlib.pyplot as plt
from layerpot import *
'''
    Test functions for module layerpot
'''
class Point:
  def __init__(self, x=(), nx=()):
    self.x = np.array(x)
    self.nx = np.array(nx)
    self.flag_inside_s = np.ones(len(x))
# -----------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------
def exactkernelK(s, index_target_point):
    '''
    This function computes exact kernel K (from double layer)
    index_target_point: it's not computed the undefined value at source point s[index_target_point]
    '''
    t = Point([s.x[index_target_point]], [s.nx[index_target_point]])
    d = t.x - s.x
    d[index_target_point] = np.nan
    r = np.abs(d)
    # ---------------------------------------------------------------
    # the two definitions can differ for floating point error
    # cosphi = np.real( np.conj(s.nx) * d) / r
    cosphi = scalar(d / r, s.nx)
    # ---------------------------------------------------------------
    # test plot: you can see instability for s = circumference
    #            when the kernel is expected constant
    if False:
        plt.figure()
        plt.title('cosphi')
        plt.plot(fundsol_deriv_negate(r = r[1:], cosphi = 1) * cosphi[1:])
        plt.plot(fundsol_deriv_negate(r = np.flip(r[1:]), cosphi = 1) * np.flip(cosphi[1:]))
    # ---------------------------------------------------------------
    return fundsol_deriv_negate(r = r, cosphi = cosphi)
def exactkernelKp(s, index_target_point):
    '''
    This function computes exact kernel Kp (from derivative of single layer)
    index_target_point: it's not computed the undefined value at source point s[index_target_point]
    '''
    t = Point([s.x[index_target_point]], [s.nx[index_target_point]])
    d = t.x - s.x
    d[index_target_point] = np.nan
    r = np.abs(d)
    # ---------------------------------------------------------------
    # the two definitions can differ for floating point error
    # cosphi = np.real( np.conj(s.nx) * d) / r
    cosphi = scalar(d / r, t.nx)
    return - fundsol_deriv_negate(r = r, cosphi = cosphi)
def exactkernelS():
    return
# -----------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------
def plot_kernel_SegmentPoint(s, index_target_point = 0, layertype = 'K'):
    '''
        This function takes a segment 's' object and the index 'k' of target point 'p' extract from s[]
        Then the approximated kernel defined for each point in s (s=s, t=p) is visualized
        - layertype: it specifies the type of kernel/layerpotential
    '''
    switcher = {'S' : layerpotS, 'K' : layerpotD, 'Kp' : layerpotSD}
    A = switcher[layertype](s=s, noweights=True) # target it's not specified to trigger the self approximation
    # --- plot ---    
    plt.figure()
    kernelfunc_approx = A[index_target_point, :]
    plt.plot(kernelfunc_approx, '1-')
    plt.plot(index_target_point, - s.kappa[index_target_point] / 4 / pi, 'r*')

    switcher_exact = {'S' : exactkernelS, 'K' : exactkernelK, 'Kp' : exactkernelKp}
    kernelfunc = switcher_exact[layertype](s, index_target_point)
    plt.figure(1)
    plt.plot(kernelfunc, '2-')