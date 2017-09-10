from src import *
import dill

# dill.load_session('runs/dill/dill_two_ellipse_delta1e-6.pkl')
# dill.load_session('runs/dill/dill_three_ellipse_delta1e-6.pkl')
dill.load_session('runs/dill/dill_sharp_ellipse_delta1e-6.pkl')
p.ipb_opt_append(4)

# dill.dump_session('runs/dill/dill_two_ellipse_delta1e-6_app.pkl')
# dill.dump_session('runs/dill/dill_three_ellipse_delta1e-6_app.pkl')
dill.dump_session('runs/dill/dill_sharp_ellipse_delta1e-6_app.pkl')
