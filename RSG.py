from pyxtal.crystal import random_crystal
from pymatgen.io.cif import CifWriter
import numpy as np
np.random.seed(seed=0)
composition_speciess = ['Ba', 'Rh', 'O']
candidates = [[1, 1, 3]]

space_number = np.random.randint(1, 230)
composition_ratio = candidates[int(np.random.choice(
    [0], 1, p=[1.0]))]
volume = np.random.choice([0.8, 1.0, 1.2], 1, p=[0.1, 0.8, 0.1])

print(space_number, composition_ratio, volume)

rand_crystal = random_crystal(
    space_number, composition_speciess, composition_ratio, volume)
try:
    CifWriter(rand_crystal.struct).write_file("hoge.cif")
except AttributeError:
    pass
