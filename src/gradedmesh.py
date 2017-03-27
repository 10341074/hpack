from numpy import pi
import numpy as np

p = 2.0

def v(s):
  return (1.0 / p - 0.5) * ( (pi - s) / pi ) ** 3 + 1.0 / p * (s - pi) / pi + 0.5

def w(s):
  return 2.0 * pi * v(s)**p / ( v(s)**p + v(2.0 * pi - s)**p )

def vp(s):
  return 3.0 * (1.0 / p - 0.5) * ( (pi - s) / pi ) ** 2 * (-1.0 / pi) + 1.0 / p / pi

def wp(s):
  v1 = v(s)
  v2 = v(2.0 * pi -s)
  vp1 = vp(s)
  vp2 = vp(2.0 * pi -s)
  den = v1 ** p + v2 ** p
  return 2.0 * pi / den**2 * (p * v1**(p-1) * vp1 * den - v1 ** p * (p * v1 ** (p-1) * vp1 - p * v2 ** (p-1) * vp2 ))

# w(s) = (2*pi*((s - pi)/(p*pi) - ((1/p - 1/2)*(s - pi)^3)/pi^3 + 1/2)^p)/((((1/p - 1/2)*(s - pi)^3)/pi^3 - (s - pi)/(p*pi) + 1/2)^p + ((s - pi)/(p*pi) - ((1/p - 1/2)*(s - pi)^3)/pi^3 + 1/2)^p)

# w_p(s) = (2*pi*(p*(1/(p*pi) - (3*(1/p - 1/2)*(s - pi)^2)/pi^3)*(((1/p - 1/2)*(s - pi)^3)/pi^3 - (s - pi)/(p*pi) + 1/2)^(p - 1) - p*(1/(p*pi) - (3*(1/p - 1/2)*(s - pi)^2)/pi^3)*((s - pi)/(p*pi) - ((1/p - 1/2)*(s - pi)^3)/pi^3 + 1/2)^(p - 1))*((s - pi)/(p*pi) - ((1/p - 1/2)*(s - pi)^3)/pi^3 + 1/2)^p)/((((1/p - 1/2)*(s - pi)^3)/pi^3 - (s - pi)/(p*pi) + 1/2)^p + (- ((1/p - 1/2)*(s - pi)^3)/pi^3 + (s - pi)/(p*pi) + 1/2)^p)^2 + (2*p*pi*(1/(p*pi) - (3*(1/p - 1/2)*(s - pi)^2)/pi^3)*((s - pi)/(p*pi) - ((1/p - 1/2)*(s - pi)^3)/pi^3 + 1/2)^(p - 1))/((((1/p - 1/2)*(s - pi)^3)/pi^3 - (s - pi)/(p*pi) + 1/2)^p + ((s - pi)/(p*pi) - ((1/p - 1/2)*(s - pi)^3)/pi^3 + 1/2)^p)
