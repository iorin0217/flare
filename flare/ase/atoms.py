import os
import sys
import inspect
from time import time
from copy import deepcopy

import numpy as np
from ase import Atoms
from flare.utils.learner import get_max_cutoff


class FLARE_Atoms(Atoms):
    """
    The `FLARE_Atoms` class is a child class of ASE `Atoms`, 
    which has completely the same usage as the primitive ASE `Atoms`, and
    in the meanwhile mimic `Structure` class. It is used in the `OTF` module
    with ASE engine (by `OTF_ASE` module). It enables attributes to be 
    obtained by both the name from ASE `Atoms` and `Structure`.
    The input arguments are the same as ASE `Atoms`.
    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.prev_positions = np.zeros_like(self.positions)

    @staticmethod
    def from_ase_atoms(atoms):
        """
        Args:
            atoms (ASE Atoms): the ase atoms to build from
        """
        kw = [
            "symbols",  # use either "symbols" or "number"
            "positions",
            "tags",
            "momenta",
            "masses",
            "magmoms",
            "charges",
            "scaled_positions",
            "cell",
            "pbc",
            "celldisp",
            "constraint",
            "calculator",
            "info",
            "velocities",
        ]
        # The keywords are either the attr of atoms, or in dict atoms.arrays
        kwargs = {}
        for key in kw:
            try:
                kwargs[key] = getattr(atoms, key)
            except:
                if key in atoms.arrays:
                    kwargs[key] = atoms.arrays[key]
        kwargs["calculator"] = atoms.calc
        return FLARE_Atoms(**kwargs)

    @property
    def nat(self):
        return len(self)

    @property
    def species_labels(self):
        return self.symbols

    @property
    def coded_species(self):
        return self.numbers

    @property
    def forces(self):
        return self.get_forces()

    @forces.setter
    def forces(self, forces_array):
        pass

    @property
    def potential_energy(self):
        return self.get_potential_energy()

    @potential_energy.setter
    def potential_energy(self, potential_energy_array):
        pass

    @property
    def stress(self):
        return self.get_stress()

    @stress.setter
    def stress(self, stress_array):
        pass

    @property
    def local_energies(self):
        try:  # when self.calc is not FLARE, there's no get_local_energies()
            local_energies = self.get_local_energies()
        except:
            local_energies = np.zeros(self.nat)
        return local_energies

    @local_energies.setter
    def local_energies(self, local_energies_array):
        pass

    @property
    def partial_stresses(self):
        try:  # when self.calc is not FLARE, there's no get_spartial_stresses()
            partial_stresses = self.get_partial_stresses()
        except:
            partial_stresses = np.zeros((self.nat, 6))
        return partial_stresses

    @partial_stresses.setter
    def partial_stresses(self, partial_stresses_array):
        pass

    @property
    def energy_stds(self):
        try:  # when self.calc is not FLARE, there's no get_energy_stds()
            energy_stds = self.get_energy_stds()
        except:
            energy_stds = 0
        return energy_stds

    @energy_stds.setter
    def energy_stds(self, energy_stds_array):
        pass

    @property
    def stress_stds(self):
        try:  # when self.calc is not FLARE, there's no get_stress_stds()
            stress_stds = self.get_stress_stds()
        except:
            stress_stds = np.zeros(6)
        return stress_stds

    @stress_stds.setter
    def stress_stds(self, stress_stds_array):
        pass

    @property
    def partial_stress_stds(self):
        try:  # when self.calc is not FLARE, there's no get_partial_stress_stds()
            partial_stress_stds = self.get_partial_stress_stds()
        except:
            partial_stress_stds = np.zeros((self.nat, 6))
        return partial_stress_stds

    @partial_stress_stds.setter
    def partial_stress_stds(self, partial_stress_stds_array):
        pass

    @property
    def local_energy_stds(self):
        try:  # when self.calc is not FLARE, there's no get_local_energy_stds()
            local_energy_stds = self.get_local_energy_stds()
        except:
            local_energy_stds = np.zeros(self.nat)
        return local_energy_stds

    @local_energy_stds.setter
    def local_energy_stds(self, local_energy_stds_array):
        pass

    @property
    def stds(self):
        try:  # when self.calc is not FLARE, there's no get_stds()
            stds = self.get_stds()
        except:
            stds = np.zeros_like(self.positions)
        return stds

    @stds.setter
    def stds(self, stds_array):
        pass

    def wrap_positions(self):
        return self.get_positions(wrap=True)

    @property
    def wrapped_positions(self):
        return self.get_positions(wrap=True)

    @property
    def max_cutoff(self):
        return get_max_cutoff(self.cell)
