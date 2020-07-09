'''
python MLE.py gp_{num}.pickle, log.txt
'''
import os
import sys
import numpy as np
import pandas as pd
from flare.gp import GaussianProcess
from flare.struc import Structure

gp_pickle = sys.argv[1]
log_txt = sys.argv[2]
gp_model = GaussianProcess.from_file(gp_pickle)
with open(log_txt, "r") as f:
    tmp = f.readlines()
    logs = [i.split("\n")[0] for i in tmp]
previous_type, previous_num = logs[-1].split(" ")

# check previous calc
if previous_type == "DFT":
    exp_path = os.path.dirname(log_txt)
    with open(f"{exp_path}/dft_targets_{previous_num}.txt", "r") as f:
        tmp = f.readlines()
        targets = [i.split("\n")[0] for i in tmp]
    # update gp db
    for i in range(0, len(targets), 2):
        structure_pickle = targets[i]
        add_targets = [int(j) for j in targets[i + 1].split(",")]
        structure_saved = pd.read_pickle(structure_pickle)
        structure = Structure(cell=structure_saved._cell,
                              positions=structure_saved._positions, species=structure_saved.species_labels)
        gp_model.update_db(structure, structure_saved.forces,
                           custom_range=add_targets, energy=structure_saved.potential_energy)
    # update kernel matrix
    gp_model.set_L_alpha()
# update hyps
if previous_num % 10 == 0:
    gp_model.train()
# save
gp_model.write_model(f'gp_{previous_num+1}', format='pickle')
print(*(logs + [f"update {previous_num}"]),
      sep="\n", end="\n", file=open(log_txt, "w"))
