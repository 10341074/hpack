from src import *
import dill

dill.load_session('runs/dill/dill_two_ellipse_delta1e-6.pkl')
p.ipb_opt_append(2)

dill.dump_session('runs/dill/dill_two_ellipse_delta1e-6_app.pkl')
