import sys
import time
from importlib import reload

import numpy as np
import matplotlib.pyplot as plt
import scipy.linalg as linalg
import numpy.linalg
#-------------------------------
import layerpot as ly
import segment as sg
import shapes as sh
import geometries as gm
import plot
import directproblem as dpb
import inverseproblem as ipb
import mainpb as  m
# import problem as pb
import linfunc as linf
import lintype as lint

import basis as bs

import tools as tls
import setups

#-------------------------------
import importlib
import types

def reloadall(module):
    importlib.reload(module)
    for child in module:
        if isinstance(child, types.ModuleType):
            reloadall(child)
