def v1(z):
  return z.real**3 - 3 * z.real * z.imag**2
def v1_p(z):
  return (3 * z.real**2 - 3 * z.imag**2) + 1j * ( - 6 * z.real * z.imag)