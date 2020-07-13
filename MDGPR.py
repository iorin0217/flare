'''
python MDGPR.py path/to/{composition0}_{num}.pickle path/to/{composition0}_{method}_{num}.pickle gp_{num}.pickle
'''
import os
import sys
import datetime
import pandas as pd
from flare.ase.atoms import FLARE_Atoms
from flare.ase.calculator import FLARE_Calculator
from flare.gp import GaussianProcess
from ase import units
from ase.md.npt import NPT

if __name__ == "__main__":
    # input
    structure_path = sys.argv[1]
    structure = pd.read_pickle(structure_path)
    outdir = os.path.dirname(structure_path)
    step_num = int(os.path.basename(
        structure_path).split(".")[0].split("_")[-1])
    comp = os.path.basename(structure_path).split("_")[0]
    compdir = os.path.dirname(outdir)
    # trust structure.cell is upper triangular matrix
    efs_path = sys.argv[2]
    efs = pd.read_pickle(efs_path)
    gp_pickle = sys.argv[3]
    gp_model = GaussianProcess.from_file(gp_pickle)
    # setup flare calc for MD
    flare_calculator_MD = FLARE_Calculator(gp_model,
                                           par=True,
                                           mgp_model=None,
                                           use_mapping=False)
    flare_calculator_MD.results = efs
    # setup MD
    # 5fs*2 -> 10ps (1000step_num)
    # 300~1500K=0.02585~0.12926eV 0~10GPa=0~0.06ev/A^3
    # sequence (300,0) -> (300,10) -> (1500,10) -> (1500,0)
    if step_num < 250:
        temperature, pressure = 0.02585, 0.0
    elif step_num < 500:
        temperature, pressure = 0.02585, 0.06
    elif step_num < 750:
        temperature, pressure = 0.12926, 0.06
    else:
        temperature, pressure = 0.12926, 0.0
    structure.set_calculator(flare_calculator_MD)
    md = NPT(atoms=structure, timestep=1 * units.fs, temperature=temperature,
             externalstress=pressure, ttime=None, pfactor=3375)
    # run MD
    md_start_time = datetime.datetime.now()
    md.run(1)
    md_end_time = datetime.datetime.now()
    # record_state
    md_log = [f"md_start_time : {md_start_time}", f"md_end_time : {md_end_time}", f"gibbs_free_energy : {md.get_gibbs_free_energy()}", f"total_energy : {structure.get_total_energy()}",
              f"kinetic_energy : {structure.get_kinetic_energy()}", f"temperature : {structure.get_temperature()}", structure.get_velocities()]
    print(*md_log, sep="\n", end="\n",
          file=open(f"{outdir}/md_log_{step_num}.txt", "w"))
    del flare_calculator_MD  # for memory save
    # run GPR for new structure
    flare_calculator_GPR = FLARE_Calculator(gp_model,
                                            par=True,
                                            mgp_model=None,
                                            use_mapping=False)
    structure.set_calculator(flare_calculator_GPR)
    structure.calc.calculate(structure)
    # save structure and efs
    new_dir = compdir + "/" + comp + "_" + str(step_num + 1)
    os.mkdir(new_dir)
    pd.to_pickle(structure, f"{new_dir}/{comp}_{step_num + 1}.pickle")
    results = structure.calc.show_results()
    pd.to_pickle(results, f"{new_dir}/{comp}_gpr_{step_num + 1}.pickle")
