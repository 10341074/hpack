def e():
  exit()
def issym(A):
  return max([max(r) for r in abs(A-A.T)])
def printmat(A, rind, cind):
  for r in A[rind]:
    for c in r[cind]:
      print(c, '\t', end='')
    print()
