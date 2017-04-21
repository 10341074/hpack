from numpy import pi, sin, cos, log

# circumference
def cZ(t, c = 0, a = 1, b = 1, aff=(0, 1)):
  return aff[0] + aff[1] * complex(c + a * cos(2 * pi * t) + 1j * b * sin(2 * pi * t))
def cZp(t, c = 0, a = 1, b = 1, aff=(0, 1)):
  return aff[1] * complex(-2 * pi * a * sin(2 * pi * t) + 1j * 2 * pi * b * cos(2 * pi * t))
def cZpp(t, c = 0, a = 1, b = 1, aff=(0, 1)):
  return aff[1] * complex(-2 * pi * 2 * pi * a * cos(2 * pi * t) - 1j * 2 * pi * 2 * pi * b * sin(2 * pi * t))

# segment
def sZ(t, p = 0, q = 1, aff=(0, 1)):
  return aff[0] + aff[1] * complex(p + t * (q - p))
def sZp(t, p = 0, q = 1, aff=(0, 1)):
  return aff[1] * complex(q - p)
def sZpp(t, p = 0, q = 1):
  return aff[1] * complex(0)

# kite
def kZ(t, a = [], aff=(0, 1)):
  return aff[0] + aff[1] * complex(cos(2*pi*t) + 0.65 * cos(2*2*pi*t) - 0.65 + 1j * 1.5 * sin(2*pi*t))
def kZp(t, a = [], aff=(0, 1)):
  return aff[1] * complex(-2*pi*sin(2*pi*t) - 2 * 0.65 * 2*pi*sin(2*2*pi*t) + 1j * 1.5 * 2*pi*cos(2*pi*t))
def kZpp(t, a = [], aff=(0, 1)):
  return aff[1] * complex(-2*pi*2*pi*cos(2*pi*t) - 2 * 2 * 0.65 * 2*pi*2*pi*cos(2*2*pi*t) - 1j * 1.5 * sin(2*pi*t))

# drop
def dZ(t, a = [], aff=(0, 1)):
  return aff[0] + aff[1] * complex(2.0 * sin(2 * pi * t / 2) - 1j * sin(2 * pi * t))
def dZp(t, a = [], aff=(0, 1)):
  return aff[1] * complex(2 * pi / 2 * 2.0 * cos(2 * pi * t / 2) - 2 * pi * 1j * cos(2 * pi * t))
def dZpp(t, a = [], aff=(0, 1)):
  return aff[1] * complex(-2 * pi * 2 * pi / 4 * 2.0 * sin(2 * pi * t / 2) + 2 * pi * 2 * pi * 1j * sin(2 * pi * t))

def circle(c = 0, r = 1):
  args=(c, r, r)
  return (cZ, cZp, cZpp, args)
def ellipse(c = 0, a = 1, b = 2):
  args=(c, a, b)
  return (cZ, cZp, cZpp, args)
def line(p = 0, q = 1):
  args=(p, q)
  return (sZ, sZp, sZpp, args)
