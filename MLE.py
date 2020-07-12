'''
python MLE.py gp_{num}.pickle, log.txt
'''
import os
import sys
import datetime
import pandas as pd
from flare.gp import GaussianProcess

if __name__ == "__main__":
    start_time = datetime.datetime.now()
    gp_pickle = sys.argv[1]
    log_txt = sys.argv[2]
    gp_model = GaussianProcess.from_file(gp_pickle)
    start_hyps = gp_model.hyps
    step_num = int(gp_pickle.split(".")[0].split("_")[1]) + 1
    exp_path = os.path.dirname(log_txt)
    with open(log_txt, "r") as f:
        tmp = f.readlines()
        logs = [i.split("\n")[0] for i in tmp]
    OTF_flag = logs[-1].split(" ")[-1]

    # check previous calc
    if OTF_flag == "False":
        with open(f"{exp_path}/dft_targets_{step_num}.txt", "r") as f:
            tmp = f.readlines()
            dft_targets = [i.split("\n")[0] for i in tmp]
        # update gp db
        for i in range(0, len(dft_targets), 2):
            structure_pickle = dft_targets[i]
            add_targets = [int(j) for j in dft_targets[i + 1][1:-1].split(",")]
            structure = pd.read_pickle(structure_pickle)
            results_pickle = structure_pickle.split(
                ".")[0][:-1 * len(str(step_num))] + f"dft_{step_num}.pickle"
            results = pd.read_pickle(results_pickle)
            gp_model.update_db(structure, results["forces"],
                               custom_range=add_targets, energy=results["energy"])
        # update kernel matrix
        gp_model.set_L_alpha()
    # update hyps
    if step_num % 50 == 0:
        gp_model.train()
    # save
    end_time = datetime.datetime.now()
    end_hyps = gp_model.hyps
    gp_model.write_model(f'gp_{step_num}', format='pickle')
    print(*(logs + [f"DFT {step_num} end_time {start_time}", f"MLE {step_num} end_time {end_time} start_hyps {start_hyps} end_hyps {end_hyps}"]),
          sep="\n", end="\n", file=open(log_txt, "w"))
