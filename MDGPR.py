'''
python MDGPR.py path/to/flare_atoms_{num}.pickle path/to/flare_atoms_{method}_{num}.pickle gp_{num}.pickle
'''
import pandas as pd
from flare.ase.calculator import FLARE_Calculator
from flare.ase.atoms import FLARE_Atoms
from flare.struc import Structure
from ase import Atoms
from ase import units
from ase.md.npt import NPT

# input
structure_path = sys.argv[1]
structure = pd.read_pickle(structure_path)
prefix = structure_path.split("_")[0]
step = int(structure_path.split("/")[-1].split(".")[0].split("_")[-1])
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
# 5fs*2 -> 10ps (1000step)
# 300~1500K=0.02585~0.12926eV 0~10GPa=0~0.06ev/A^3
# sequence (300,0) -> (300,10) -> (1500,10) -> (1500,0)
if step < 250:
    temperature, pressure = 0.02585, 0.0
elif step < 500:
    temperature, pressure = 0.02585, 0.06
elif step < 750:
    temperature, pressure = 0.12926, 0.06
else:
    temperature, pressure = 0.12926, 0.0
structure.set_calculator(flare_calculator_MD)
md = NPT(atoms=structure, timestep=5 * units.fs, temperature=temperature,
         externalstress=pressure, ttime=None, pfactor=3375)
# run MD
md.run(2)
# record_state

del flare_calculator_MD  # for memory save
# run GPR for new structure
flare_calculator_GPR = FLARE_Calculator(gp_model,
                                        par=True,
                                        mgp_model=None,
                                        use_mapping=False)
structure.set_calculator(flare_calculator_GPR)
structure.calc.calculate(structure)
# save structure and efs
pd.to_pickle(structure, f"{prefix}_{step + 1}.pickle")
results = structure.calc.show_results
pd.to_pickle(results, f"{prefix}_gpr_{step + 1}.pickle")
