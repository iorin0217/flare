from scipy.sparse.linalg import LinearOperator, svds
import numpy as np
from dscribe.descriptors import SOAP
soap = SOAP(
    species=[56, 8, 45],
    sigma=0.0875,
    rcut=10.5,
    nmax=10,
    lmax=9,
    periodic=True,
    sparse=False,
    average=False
)


def descriptor_svd(at_descs, num, do_vectors='vh'):
    def mv(v):
        return np.dot(at_descs, v)

    def rmv(v):
        return np.dot(at_descs.T, v)

    A = LinearOperator(at_descs.shape, matvec=mv, rmatvec=rmv, matmat=mv)

    return svds(A, k=num, return_singular_vectors=do_vectors)


soap_B1R1O3_1 = soap.create(B1R1O3_1)
soap_vecs = []
all_dict = {}
envs = [(B1R1O3_1, soap_B1R1O3_1)]
for (crystal, soap) in [envs[1]]:
    index_list = []
    for index, element in enumerate(crystal.get_chemical_symbols()):
        if element == "Rh":
            soap_vecs.append(soap[index])
            index_list.append(index)
            all_dict[crystal.get_chemical_formula()] = index_list

at_descs = np.array(soap_vecs).T

num = 3
m = at_descs

(u, s, vt) = descriptor_svd(m, min(max(1, int(num/2)), min(m.shape)-1))
c_scores = np.sum(vt**2, axis=0)/vt.shape[0]
selected = sorted(np.argsort(c_scores)[-num:])
