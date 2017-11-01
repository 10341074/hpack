from src import *

def plot(no = 80, nb = 80):
  t = 2 * np.pi / no * np.arange(no)
  # B = bs.basis_lagr_trig(no, t)
  B = bs.basis_lagr(no, t)
  print(B.shape)
  plt.plot(B[:,:6],'*-')
  plt.show(block=False)

  p = m.EIT()
  p.domain()
  plt.figure()
  p.plot_domain()
  p.solver()

  B = p.K.dot(B)

  plt.figure()
  plt.plot(B[:,:6],'*-')
  plt.show(block=False)

  return

if __name__ == "__main__":
  plot()
  inp = input("Press a key")
