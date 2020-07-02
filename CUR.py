from quippy import descriptors, Atoms
import sys
import os
import numpy as np
from scipy.sparse.linalg import LinearOperator, svds
import ase.io


def file_to_file(in_file, out_file, method_name, num, prev_config_files=[], method_kwargs=None):
    print "selecting by descriptor from ", in_file, "to", out_file
    # read configs
    ats = ase.io.read(in_file, ":")

    # create list of previous descriptors from files in glob
    descs_prev_list = []
    for prev_file in prev_config_files:
        descs_prev_list.extend([at.info["descriptor_vec"]
                                for at in ase.io.read(prev_file, ":")])

    print "number of configs", len(ats), "previous descriptors", len(descs_prev_list)

    selected_ats = globals()[method_name](
        ats, descs_prev_list, num, **method_kwargs)

    # write with descriptor vectors, in case that's needed at later iterations
    ase.io.write(os.path.join(os.path.dirname(out_file),
                              'descriptor_vec.'+os.path.basename(out_file)), selected_ats)

    # remove descriptor and write again
    for at in selected_ats:
        del at.info["descriptor_vec"]
    ase.io.write(out_file, selected_ats)


def calc_descriptors(ats, descriptor_file=None, descriptor_string=None):
    if descriptor_file is None and descriptor_string is None:
        desc = None  # hope that it's set in atoms objects
    else:
        if descriptor_file is not None:
            if descriptor is not None:
                raise Exception(
                    'descriptor and descriptor_file both specified')
            with open(descriptor_file, "r") as f:
                descriptor_string = f.readline()
        desc = descriptors.Descriptor(descriptor_string+" average")
        desc_cutoff = desc.cutoff()

    sys.stderr.write("descriptors file "+str(descriptor_file)+" string "+str(
        descriptor_string)+" has descriptor_vec " + str("descriptor_vec" in ats[0].info)+"\n")
    for (at_i, at) in enumerate(ats):
        if "descriptor_vec" not in at.info:
            at_q = Atoms(at)
            at_q.set_cutoff(desc_cutoff)
            at_q.calc_connect()
            at.info["descriptor_vec"] = desc.calc(at_q).descriptor[0]
            if at_i % 100 == 99:
                sys.stderr.write('{}'.format(int(at_i/100) % 10))
    sys.stderr.write('\n')


def descriptor_svd(at_descs, num, do_vectors='vh'):
    def mv(v):
        return np.dot(at_descs, v)

    def rmv(v):
        return np.dot(at_descs.T, v)

    A = LinearOperator(at_descs.shape, matvec=mv, rmatvec=rmv, matmat=mv)

    return svds(A, k=num, return_singular_vectors=do_vectors)


def CUR(ats, descs_prev_list, num, stochastic=True, kernel_exp=0.0):
    # column vectors of descriptors
    at_descs = np.array([at.info["descriptor_vec"] for at in ats]).T

    # do SVD on kernel if desired
    if kernel_exp > 0.0:
        m = np.matmul(at_descs.T, at_descs)**kernel_exp
    else:
        m = at_descs

    print "len(ats)", len(ats), "at_descs.shape", at_descs.shape, "m.shape", m.shape

    (u, s, vt) = descriptor_svd(m, min(max(1, int(num/2)), min(m.shape)-1))

    c_scores = np.sum(vt**2, axis=0)/vt.shape[0]
    # sys.stderr.write("c_scores indices "+ str(np.argsort(c_scores))+"\n")
    # sys.stderr.write("c_scores values "+ str(np.sort(c_scores))+"\n")
    if stochastic:
        selected = sorted(np.random.choice(
            range(len(ats)), size=num, replace=False, p=c_scores))
    else:
        selected = sorted(np.argsort(c_scores)[-num:])

    return [ats[i] for i in selected]


def greedy_farthest_point(ats, descs_prev_list, num, kernel_exp):
    # row vectors of descriptors
    at_descs = np.array([at.info["descriptor_vec"] for at in ats])
    print "at_descs.shape", at_descs.shape

    cur_similarities = np.matmul(at_descs, at_descs.T)**kernel_exp
    # zero diagonal so max similarity ignores similarity to self
    cur_similarities -= np.diag(np.diag(cur_similarities))
    print "cur_similarities.shape", cur_similarities.shape,

    # in prev similarities, row corresponds to available point, and column to all past points
    if len(descs_prev_list) > 0:
        # column vectors
        descs_prev = np.array(descs_prev_list).T
        prev_similarities = np.matmul(at_descs, descs_prev)**kernel_exp
        print "prev_similarities.shape", prev_similarities.shape,
    else:
        prev_similarities = np.zeros((len(ats), 1))
    print ""

    selected_similarities = np.zeros((len(ats), num))
    selected = []
    for i in range(num):
        if i == 0 and len(descs_prev_list) == 0:
            # if nothing previous, first point is random
            selected_ind = np.random.randint(cur_similarities.shape[0])
            print "random selected_ind", selected_ind, "out of", cur_similarities.shape[0]
        else:
            # for each available point, first find similarity to previous point that is closest to it (highest similarity)
            # then pick available point that is most different (lowest similarity among points) from whatever it's most similar to
            # print "max similarity among selected",np.max(selected_similarities, axis=1)
            # print "max similarity among prev",np.max(prev_similarities, axis=1)
            max_similarity = np.max(np.stack((np.max(selected_similarities, axis=1), np.max(
                prev_similarities, axis=1)), axis=1), axis=1)
            # print "max similarity overall",max_similarity
            selected_ind = np.argmin(max_similarity)
            print "minmax selected_ind", selected_ind, "dist", max_similarity[selected_ind]

        # select this point and copy its similarity to each other current point to selected_similarities
        selected.append(selected_ind)
        selected_similarities[:, i] = cur_similarities[:, selected_ind]
        # set this point's similarity to all other points (this point's row) to be very high, so max_similarity ends up high, and argmin never picks it
        prev_similarities[selected_ind, :] = np.finfo('float64').max
        selected_similarities[selected_ind, :] = np.finfo('float64').max

    return [ats[i] for i in selected]


def traj_MDS_distance(ats, MDS_basis, spacing, include_0=False, metric=None):
    if metric is None:
        metric = np.identity(MDS_basis.shape[0])

    selected_ats = []
    if include_0:
        selected_ats.append(ats[0])

    at_prev = ats[0]
    p_prev = np.dot(at_prev.info['descriptor_vec'], MDS_basis.T)
    for at in ats[1:]:
        p = np.dot(at.info['descriptor_vec'], MDS_basis.T)
        if np.dot((p-p_prev), np.dot(metric, (p-p_prev))) > spacing**2:
            selected_ats.append(at)
            at_prev = at
            p_prev = p

    return selected_ats
