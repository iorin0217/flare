'''
python DFT.py path/to/{composition0}_{num}.pickle
'''
import pandas as pd
import os
import sys
from flare.ase.atoms import FLARE_Atoms
from ase.calculators.espresso import Espresso

if __name__ == "__main__":
    # input
    input_path = sys.argv[1]
    label = input_path.split("/")[-2]
    input_dir = os.path.dirname(input_path)
    pseudo_dir = os.path.dirname(os.path.dirname(input_dir))
    structure = pd.read_pickle(input_path)
    # os config
    no_cpus = 32
    npool = 4
    input_file_name = input_dir + "/" + label + ".pwi"
    output_file_name = input_dir + "/" + label + ".pwo"
    pw_loc = os.environ.get("PWSCF_COMMAND")
    os.environ["ASE_ESPRESSO_COMMAND"] = f"mpirun -np {no_cpus} {pw_loc} -npool {npool} < {input_file_name} > {output_file_name}"
    # setup calculator
    threshold = len(structure) * 0.2e-9
    input_data = {"control": {"pseudo_dir": f"{pseudo_dir}/", "calculation": "scf"},
                  "system": {"ibrav": 0, "ecutwfc": 50, "ecutrho": 400, "occupations": "tetrahedra_opt"}, "electrons": {"conv_thr": threshold, "mixing_beta": 0.2}}
    pseudo_files = {"Ba": "Ba.pbe-spn-kjpaw_psl.1.0.0.UPF",
                    "Rh": "Rh.pbe-spn-kjpaw_psl.1.0.0.UPF", "O": "O.pbe-n-kjpaw_psl.0.1.UPF"}
    dft_calc = Espresso(pseudopotentials=pseudo_files, label=input_dir + "/" + label,
                        tstress=True, tprnfor=True, nosym=True, input_data=input_data, kpts=(6, 6, 6))
    structure.set_calculator(dft_calc)
    # get results
    results = {}
    results["energy"] = structure.get_potential_energy()
    results["forces"] = structure.get_forces()
    results["stress"] = structure.get_stress()
    # save
    comp, num = label.split("_")
    pd.to_pickle(results, f"{input_dir}/{comp}_dft_{num}.pickle")
