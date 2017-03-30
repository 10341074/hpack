# execfile('../hpack/src/__load__.py')
import matplotlib.pyplot as plt
import numpy as np
import scipy.linalg as lin

# python3
from importlib import reload

# __modules__ are scripts, import *
# modules contains classes, functions, __main__ to be tested
import __types__
import __data__
import __tools__ # for help(__tools__)
from __tools__ import *

import layerpot as ly
import segment as sg
import plot
import directproblem as dpb
import inverseproblem as ipb
