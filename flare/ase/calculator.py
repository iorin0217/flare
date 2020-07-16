""":class:`FLARE_Calculator` is a calculator compatible with `ASE`.
You can build up `ASE Atoms` for your atomic structure, and use `get_forces`,
`get_potential_energy` as general `ASE Calculators`, and use it in
`ASE Molecular Dynamics` and our `ASE OTF` training module. For the usage
users can refer to `ASE Calculator module <https://wiki.fysik.dtu.dk/ase/ase/calculators/calculators.html>`_
and `ASE Calculator tutorial <https://wiki.fysik.dtu.dk/ase/ase/atoms.html#adding-a-calculator>`_."""

import warnings
import numpy as np
import multiprocessing as mp
from flare.env import AtomicEnvironment
from flare.struc import Structure
from flare.mgp import MappedGaussianProcess
from flare.predict import (
    predict_on_structure_par_en,
    predict_on_structure_en,
    predict_on_structure_efs,
    predict_on_structure_efs_par,
)
from ase.calculators.calculator import Calculator


class FLARE_Calculator(Calculator):
    """
    Build FLARE as an ASE Calculator, which is compatible with ASE Atoms and
    Molecular Dynamics.
    Args:
        gp_model (GaussianProcess): FLARE's Gaussian process object
        mgp_model (MappedGaussianProcess): FLARE's Mapped Gaussian Process
            object. `None` by default. MGP will only be used if `use_mapping`
            is set to True.
        par (Bool): set to `True` if parallelize the prediction. `False` by
            default.
        use_mapping (Bool): set to `True` if use MGP for prediction. `False`
            by default.
    """

    def __init__(self, gp_model, mgp_model=None, par=False, use_mapping=False):
        super().__init__()  # all set to default values, TODO: change
        self.mgp_model = mgp_model
        self.gp_model = gp_model
        self.use_mapping = use_mapping
        self.par = par
        self.results = {}

    def get_property(self, name, atoms=None, allow_calculation=True):
        if name not in self.results.keys():
            if not allow_calculation:
                return None
            self.calculate(atoms)
        return self.results[name]

    def get_potential_energy(self, atoms=None, force_consistent=False):
        return self.get_property("energy", atoms)

    def get_forces(self, atoms):
        return self.get_property("forces", atoms)

    def get_stress(self, atoms):
        return self.get_property("stress", atoms)

    def get_energy_stds(self, atoms):
        return self.get_property("energy_stds", atoms)

    def get_stds(self, atoms):
        return self.get_property("stds", atoms)

    def get_stress_stds(self, atoms):
        return self.get_property("stress_stds", atoms)

    def get_local_energies(self, atoms):
        return self.get_property("local_energies", atoms)

    def get_partial_stresses(self, atoms):
        return self.get_property("partial_stresses", atoms)

    def get_local_energy_stds(self, atoms):
        return self.get_property("local_energy_stds", atoms)

    def get_partial_stress_stds(self, atoms):
        return self.get_property("partial_stress_stds", atoms)

    def calculate(self, atoms):
        """
        Calculate properties including: energy, local energies, forces,
            stress, uncertainties.
        Args:
            atoms (FLARE_Atoms): FLARE_Atoms object
        """

        if self.use_mapping:
            self.calculate_mgp(atoms)
        else:
            self.calculate_gp(atoms)
    '''
    def calculate_gp(self, atoms):
        # Compute energy, forces, and stresses and their uncertainties
        if self.par:
            res = predict_on_structure_efs_par(
                atoms, self.gp_model, write_to_structure=False
            )
        else:
            res = predict_on_structure_efs(
                atoms, self.gp_model, write_to_structure=False
            )

        # Set the energy, force, and stress attributes of the calculator.
        res_name = [
            "local_energies",
            "forces",
            "partial_stresses",
            "local_energy_stds",
            "force_stds",
            "partial_stress_stds",
        ]
        res_dims = [1, 3, 6, 1, 3, 6]
        for i in range(len(res_name)):
            # assert (res[i].shape[1] == res_dims[i], "shape doesn't match")
            self.results[res_name[i]] = res[i]
    '''

    def calculate_gp(self, atoms):
        if self.par:
            res = predict_on_structure_efs_par(
                atoms, self.gp_model, n_cpus=32, write_to_structure=False)
        else:
            res = predict_on_structure_efs(
                atoms, self.gp_model, write_to_structure=False)
        current_volume = np.linalg.det(np.array(atoms.cell))
        flare_stress = np.sum(res[2], 0) / current_volume
        flare_stress_stds = \
            (np.sqrt(np.sum(res[5]**2, 0)) / current_volume)
        self.results["energy"] = np.sum(res[0])
        self.results["forces"] = res[1]
        self.results["stress"] = -np.array([flare_stress[0], flare_stress[3], flare_stress[5],
                                            flare_stress[4], flare_stress[2], flare_stress[1]])  # ASE format
        self.results["energy_stds"] = np.sqrt(np.sum(res[3]**2))
        self.results["stds"] = res[4]
        self.results["stress_stds"] = np.array([flare_stress_stds[0], flare_stress_stds[3], flare_stress_stds[5],
                                                flare_stress_stds[4], flare_stress_stds[2], flare_stress_stds[1]])
        self.results["local_energies"] = res[0]
        self.results["partial_stresses"] = res[2]
        self.results["local_energy_stds"] = res[3]
        self.results["partial_stress_stds"] = res[5]

    def show_results(self):
        return self.results

    def calculate_mgp(self, atoms):
        nat = len(atoms)

        self.results["forces"] = np.zeros((nat, 3))
        self.results["partial_stresses"] = np.zeros((nat, 6))
        self.results["stds"] = np.zeros((nat, 3))
        self.results["local_energies"] = np.zeros(nat)

        rebuild_dict = {}
        newbound_dict = {}
        repredict_atoms = []
        for n in range(nat):
            chemenv = AtomicEnvironment(
                atoms, n, self.gp_model.cutoffs, cutoffs_mask=self.mgp_model.hyps_mask
            )

            # TODO: Check that stress is being calculated correctly.
            try:
                f, v, vir, e = self.mgp_model.predict(chemenv)
                self.results["forces"][n] = f
                self.results["partial_stresses"][n] = vir
                self.results["stds"][n] = np.sqrt(np.absolute(v))
                self.results["local_energies"][n] = e

            except ValueError as err_msg:  # if lower_bound error is raised
                warnings.warn("Re-build map with a new lower bound")
                re_dict = err_msg.args[0]
                nb_dict = err_msg.args[1]
                for xb in re_dict:  # collect two & three body maps
                    if xb in rebuild_dict:
                        # collect all species
                        for s_ind, spc in enumerate(re_dict[xb]):
                            if spc in rebuild_dict[xb]:
                                spc_ind = rebuild_dict[xb].index(spc)
                                if nb_dict[xb][s_ind] < newbound_dict[xb][spc_ind]:
                                    newbound_dict[xb][spc_ind] = nb_dict[xb][s_ind]
                            else:
                                rebuild_dict[xb].append(spc)
                                newbound_dict[xb].append(nb_dict[xb][s_ind])
                    else:
                        rebuild_dict[xb] = re_dict[xb]
                        newbound_dict[xb] = nb_dict[xb]

                repredict_atoms.append((n, chemenv))

        if len(rebuild_dict) > 0:
            # rebuild map for those problematic species
            for xb in rebuild_dict:
                for s_ind, spc in enumerate(rebuild_dict[xb]):
                    map_ind = self.mgp_model.maps[xb].find_map_index(spc)
                    rebuild_map = self.mgp_model.maps[xb].maps[map_ind]
                    rebuild_map.set_bounds(
                        newbound_dict[xb][s_ind], rebuild_map.bounds[1]
                    )
                    rebuild_map.build_map_container()
                    rebuild_map.build_map(self.gp_model)

            # re-predict forces, energies, etc. for those problematic atoms
            for ra in repredict_atoms:
                n, chemenv = ra
                f, v, vir, e = self.mgp_model.predict(chemenv)
                self.results["forces"][n] = f
                self.results["partial_stresses"][n] = vir
                self.results["stds"][n] = np.sqrt(np.absolute(v))
                self.results["local_energies"][n] = e

        # get global properties
        volume = atoms.get_volume()
        total_stress = np.sum(self.results["partial_stresses"], axis=0)
        self.results["stress"] = total_stress / volume
        self.results["energy"] = np.sum(self.results["local_energies"])

    def calculation_required(self, atoms, quantities):
        return True
