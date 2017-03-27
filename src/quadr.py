import numpy as np
from numpy.fft import ifft
from numpy import pi

def kress_Rjn(n):
  # function Rjn = kress_Rjn(n)
  #
  # return 2n length vector of R_j^(n) for j=0...2n-1. Takes O(n ln n) work
  # and O(n) storage, using fft for the trig sum.
  # Note the full R_{|i-j|}^(n) matrix is then circulant(kress_Rjn(N/2)).
  # See Kress MCM 1991 paper or Lin Int Eqn book p. 210
  # barnett 2/6/08

  #  Copyright (C) 2008, 2009, Alex Barnett, Timo Betcke

  if (2 * n) % 2 == 1:
    print('kress_Rjn: N=2n must be even!')
  m = np.arange(1, n)
  Rjn = -2. * pi * ifft(np.concatenate(( [0], 1.0/m, [1.0/n], 1.0/m[::-1] )))
  return np.real(Rjn)
