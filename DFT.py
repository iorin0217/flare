'''
python DFT.py path/to/ase_structure.pickle
'''
import pandas as pd
import os
import sys
from ase.calculators.espresso import Espresso
from flare.struc import Structure

# input
input_path = sys.argv[1]
tmp = input_path.split("/")
label = tmp[-2]
input_dir = input_path[:-1 * len(label + "_ase.pickle")]
pseudo_dir = "/" + tmp[1] + "/" + tmp[2] + "/" + tmp[3] + "/" + tmp[4] + "/"
structure_ase = pd.read_pickle(input_path)
# os config
no_cpus = 32
npool = 8
input_file_name = input_dir + label + ".pwi"
output_file_name = input_dir + label + ".pwo"
pw_loc = os.environ.get("PWSCF_COMMAND")
os.environ["ASE_ESPRESSO_COMMAND"] = f"mpirun -np {no_cpus} {pw_loc} -npool {npool} < {input_file_name} > {output_file_name}"
# setup calculator
input_data = {"control": {"pseudo_dir": f"{pseudo_dir}", "calculation": "scf"},
              "system": {"ibrav": 0, "ecutwfc": 50, "ecutrho": 400, "occupations": "tetrahedra_opt"}, "electrons": {"conv_thr": 1.0e-12, "mixing_beta": 0.2}}
pseudo_files = {"Ba": "Ba.pbe-spn-kjpaw_psl.1.0.0.UPF",
                "Rh": "Rh.pbe-spn-kjpaw_psl.1.0.0.UPF", "O": "O.pbe-n-kjpaw_psl.0.1.UPF"}
dft_calc = Espresso(pseudopotentials=pseudo_files, label=input_dir + label,
                    tstress=True, tprnfor=True, nosym=True, input_data=input_data, kpts=(8, 8, 8))
structure_ase.set_calculator(dft_calc)
# ASE to FLARE structure
structure_flare = Structure.from_ase_atoms(structure_ase)
structure_flare.potential_energy = structure_ase.get_potential_energy()
structure_flare.forces = structure_ase.get_forces()
structure_flare.stress = structure_ase.get_stress()
# save
pd.to_pickle(structure_flare, input_dir + label + "_flare.pickle")
