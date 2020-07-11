'''
python OTF.py md_targets_{num}.txt, log.txt
'''
import os
import sys
import time
import pandas as pd

# config
std_tolerance = 2
noise =
max_atoms_added = 12
# input
md_targets_txt = sys.argv[1]
step_num = int(os.path.basename(
    md_targets_txt).split(".")[0].split("_")[-1]) + 1
log_txt = sys.argv[2]
expdir = os.path.dirname(md_targets_txt)
# collelct
with open(md_targets_txt, "r") as f:
    tmp = f.readlines()
    md_targets = [i.split("\n")[0] for i in tmp]
# read step + 1 target
structures = []
resultss = []
for previous_structure_path in md_targets:
    compdir = os.path.dirname(os.path.dirname(previous_structure_path))
    comp = os.path.basename(compdir)
    outdir = comdpir + "/" + comp + "_" + step_num
    structure = pd.read_pickle(
        outdir + "/" + comp + "_" + step_num + ".pickle")
    results = pd.read_pickle(outdir + "/" + comp + "_" +
                             "gpr" + "_" + step_num + ".pickle")
    structures.append(structure)
    resultss.append(results)
flag, target_atoms = is_std_in_bound_par(
    std_tolerance, noise, structures, resultss, max_atoms_added)
if flag:
    # log
    # MDGPRもかいてあげる
    # dft_targets空
else:
    dft_targets.txt
    results変える

# target_atomsはtarge_structureとセット


def is_std_in_bound_par(std_tolerance, noise, structures, resultss, max_atoms_added):
    """
    Given an uncertainty tolerance and a structure decorated with atoms,
    species, and associated uncertainties, return those which are above a
    given threshold, agnostic to species.

    If std_tolerance is negative, then the threshold used is the absolute
    value of std_tolerance.

    If std_tolerance is positive, then the threshold used is
    std_tolerance * noise.

    If std_tolerance is 0, then do not check.

    :param std_tolerance: If positive, multiply by noise to get cutoff. If
        negative, use absolute value of std_tolerance as cutoff.
    :param noise: Noise variance parameter
    :param structure: Input structure
    :type structure: FLARE Structure
    :param max_atoms_added: Maximum # of atoms to add
    :return: (True,[-1]) if no atoms are above cutoff, (False,[...]) of the
            top `max_atoms_added` uncertainties
    """
    # set uncertainty threshold
    if std_tolerance == 0:
        return True, [-1]
    elif std_tolerance > 0:
        threshold = std_tolerance * np.abs(noise)
    else:
        threshold = np.abs(std_tolerance)

    # sort max stds
    nat = len(structure)
    max_stds = np.zeros((nat))
    for atom, std in enumerate(structure.stds):
        max_stds[atom] = np.max(std)
    stds_sorted = np.argsort(max_stds)
    target_atoms = list(stds_sorted[-max_atoms_added:])

    # if above threshold, return atom
    if max_stds[stds_sorted[-1]] > threshold:
        return False, target_atoms
    else:
        return True, [-1]
