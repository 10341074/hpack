from src import *
import dill
# ----------------------------------------
'''
filename is specified without format .pkl
''' 
filename = 'prova'
# ----------------------------------------
dill.load_session(filename + '.pkl')
p.ipb_opt(2)
dill.dump_session(filename + '_out.pkl')
