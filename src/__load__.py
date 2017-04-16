# execfile('../hpack/src/__load__.py')
import matplotlib.pyplot as plt
import numpy as np
import scipy.linalg as linalg
import numpy.linalg
import time

# python3
from importlib import reload

# __modules__ are scripts, import *
# modules contains classes, functions, __main__ to be tested
import __tools__ # for help(__tools__)
from __tools__ import *
from __types__ import *

import layerpot as ly
import segment as sg
import shapes as sh
import plot
import directproblem as dpb
import inverseproblem as ipb

import linfunc as linf
import lintype as lint
