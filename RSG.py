# https://pyxtal.readthedocs.io/en/latest/Usage.html
from pyxtal.crystal import random_crystal
my_crystal = random_crystal(99, ['Ba', 'Ti', 'O'], [1, 1, 3], 1.0)
my_crystal.struct
