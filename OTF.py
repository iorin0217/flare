'''
python OTF.py md_targets_{num-1}.txt, log.txt
'''
import os
import re
import sys
import subprocess
import datetime
import pandas as pd
import numpy as np
from dscribe.descriptors import SOAP
from scipy.sparse.linalg import LinearOperator, svds
from flare.ase.atoms import FLARE_Atoms


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
    :param structures: [FLARE ASE Atoms]
    :param resultss: [results]
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
    stdss = []
    indexes = []
    species = []
    for i, structure in enumerate(structures):
        stds = np.array([np.max(std_vec) for std_vec in resultss[i]["stds"]])
        stdss.extend(stds)
        index = [(i, j) for j in range(len(structure))]
        indexes.extend(index)
        species.extend(list(structure.numbers))
    # sort max stds
    stds_sorted = np.argsort(np.array(stdss))[::-1]
    # if above threshold, return atom
    if stdss[stds_sorted[0]] > threshold:
        candidates = list(
            filter(lambda i: stdss[i] > threshold, list(stds_sorted)))
        # CUR select
        species = list(set(species))
        cur_distribution = cur_select(species, structures)
        counter = 0
        targets_index = []
        for candidate in candidates:
            if counter > max_atoms_added:
                break
            else:
                index = indexes[candidate]
                ratio = cur_distribution[index]
                accept = np.random.choice([True, False], p=(ratio, 1 - ratio))
                if accept:
                    targets_index.append(index)
                    counter += 1
        return False, targets_index
    else:
        return True, []


def cur_select(species, structures):
    cur_distribution = {}
    soap = SOAP(species=species, sigma=0.0875, rcut=10.5,
                nmax=10, lmax=9, periodic=True, sparse=False, average=False)
    # decompose to specie type
    soap_species = {}
    soap_species_indexes = {}
    for specie in species:
        soap_species[specie] = []
        soap_species_indexes[specie] = []
    for i, structure in enumerate(structures):
        soap_structure = soap.create(structure)
        for j, atom in enumerate(structure):
            soap_species[atom.number].append(soap_structure[j])
            soap_species_indexes[atom.number].append((i, j))
    # cur select in each specie type
    for specie in species:
        at_descs = np.array(soap_species[specie]).T
        m = at_descs
        k = len(structures)
        (u, s, vt) = descriptor_svd(m, k)
        c_scores = np.sum(vt**2, axis=0) / vt.shape[0]
        # modify distribution to up accept frequency
        frequent = c_scores * (1 / np.max(c_scores))
        for n, p in enumerate(frequent):
            cur_distribution[soap_species_indexes[specie]
                             [n]] = np.clip(p, None, 1)
    return cur_distribution


def descriptor_svd(at_descs, num, do_vectors='vh'):
    # for sparse matrix
    def mv(v):
        return np.dot(at_descs, v)

    def rmv(v):
        return np.dot(at_descs.T, v)

    A = LinearOperator(at_descs.shape, matvec=mv, rmatvec=rmv, matmat=mv)

    return svds(A, k=num, return_singular_vectors=do_vectors)


if __name__ == "__main__":
    start_time = datetime.datetime.now()
    # config
    std_tolerance = 2
    max_atoms_added = 8
    # input
    md_targets_txt = sys.argv[1]
    expdir = os.path.dirname(md_targets_txt)
    step_num = int(os.path.basename(
        md_targets_txt).split(".")[0].split("_")[-1]) + 1
    log_txt = sys.argv[2]
    with open(log_txt, "r") as f:
        tmp = f.readlines()
        logs = [i.split("\n")[0] for i in tmp]
    noise = float(logs[-1][:-1].split(",")[-1])
    # collelct
    with open(md_targets_txt, "r") as f:
        tmp = f.readlines()
        md_targets = [i.split("\n")[0] for i in tmp]
    # read target
    structures = []
    resultss = []
    next_md_targets = []
    for i in range(0, len(md_targets), 2):
        previous_structure_path = md_targets[i]
        compdir = os.path.dirname(os.path.dirname(previous_structure_path))
        comp = os.path.basename(compdir)
        outdir = compdir + "/" + comp + "_" + str(step_num)
        structure_path = outdir + "/" + comp + "_" + str(step_num) + ".pickle"
        structure = pd.read_pickle(structure_path)
        results_path = outdir + "/" + comp + "_" + \
            "gpr" + "_" + str(step_num) + ".pickle"
        results = pd.read_pickle(results_path)
        structures.append(structure)
        next_md_targets.append(structure_path)
        resultss.append(results)
        next_md_targets.append(results_path)
    # check
    flag, targets_index = is_std_in_bound_par(
        std_tolerance, noise, structures, resultss, max_atoms_added)
    # switch
    dft_targets = []
    if not flag:
        structure_indexes = list(set([target_index[0]
                                      for target_index in targets_index]))
        for structure_index in structure_indexes:
            atom_indexes = [target_index[1]
                            for target_index in targets_index if target_index[0] == structure_index]
            dft_targets.append(next_md_targets[structure_index * 2])
            dft_targets.append(atom_indexes)
            next_md_targets[structure_index * 2 + 1] = re.sub(
                "_gpr_", "_dft_", next_md_targets[structure_index * 2 + 1])

    # write
    end_time = datetime.datetime.now()
    print(*dft_targets, sep="\n", end="\n",
          file=open(f"{expdir}/dft_targets_{step_num}.txt", "w"))
    print(*next_md_targets, sep="\n", end="\n",
          file=open(f"{expdir}/md_targets_{step_num}.txt", "w"))
    print(*(logs + [f"MDGPR {step_num - 1} end_time {start_time}", f"OTF {step_num} end_time {end_time} {flag}"]),
          sep="\n", end="\n", file=open(log_txt, "w"))
    # submit jobs
    calc_size = len(dft_targets) / 2
    if calc_size == 0:
        subprocess.run(["qsub", f"job_MLE_{step_num}.sh"])
    else:
        job_id = subprocess.run(
            ["qsub", "-J", f"1-{calc_size}", f"job_DFT_{step_num}.sh"], encoding='utf-8', stdout=subprocess.PIPE)
        subprocess.run(
            ["qsub", "-W", f"depend=afterok:{job_id.stdout.split('.')[0]}", f"job_MLE_{step_num}.sh"])
