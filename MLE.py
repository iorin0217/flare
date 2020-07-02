'''
python MLE.py gp{num}.pickle, datas.npy, log.txt
'''
import sys
import numpy as np
from flare.gp import GaussianProcess

gp_pickle = sys.argv[1]
step = gp_pickle.split("/")[-1][2:-7]
gp_model = GaussianProcess.from_file(gp_pickle)

datas_npy = sys.argv[2]
datas = np.load(datas_npy)

with open(sys.argv[3], "r") as f:
    log = f.readlines()

data_add =
# structure, dft force

# update gp db
# loop必要
if data_add not None:
    gp.update_db(*data_add)

    gp.set_L_alpha()
# update gp model
if step % 10 == 0:
    gp.train()
gp.write_model(f'gp{step+1}', format='pickle')

# log.txt書く
