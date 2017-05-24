def e():
  exit()
  return

def issym(A):
  return max([max(r) for r in abs(A-A.T)])

def mat_print(A, rind, cind):
  for r in A[rind]:
    for c in r[cind]:
      print(c, '\t', end='')
    print()
  return

def mat_max(A):
  return max([max(abs(r)) for r in A])
