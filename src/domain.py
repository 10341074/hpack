#!../../venv/bin/python
import matplotlib.pyplot as plt
import numpy as np
from scipy.linalg import solve

import my_type
import layerpot as ly
import segment as sg
import plot
import directproblem as dpb


h = 20
n = 20
c = 1. * (h + 1)/(h-1)

ns = n

rb  = 10
nsb = rb * ns
sb = sg.Segment(nsb, f=sg.circle, inargs=(0, rb))

ro  = 6
nso = ro * ns
so = sg.Segment(nso, f=sg.circle, inargs=(0, ro))

rd  = 1
nsd = rd * ns
sd = sg.Segment(nsd, f=sg.circle, inargs=(0, rd))
