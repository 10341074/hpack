def e():
  """Exit function"""
  exit()
  return
# -------------------------------------------
def issym(A):
  """A matrix: it returns norm of asymmetric part (to check A symmetry)"""
  return max([max(r) for r in abs(A - A.T)])

mat_max_asym = issym

def mat_print(A, rind, cind):
  """A matrix: it prints A"""
  for r in A[rind]:
    for c in r[cind]:
      print(c, '\t', end='')
    print()
  return

def mat_max(A):
  """A matrix: it returns """
  return max([max(abs(r)) for r in A])
# -------------------------------------------
