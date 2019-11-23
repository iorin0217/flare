"""
The cutoffs module gives a few different options for smoothly sending the GP
kernel to zero near the boundary of the cutoff sphere.
"""
from math import cos, sin, pi
from numba import njit


@njit
def hard_cutoff(r_cut, ri, ci):
    """A hard cutoff that assigns a value of 1 to all interatomic distances.
    
    :param r_cut: Cutoff value (in angstrom).
    :type r_cut: float
    :param ri: Interatomic distance.
    :type ri: float
    :param ci: Cartesian coordinate divided by the distance.
    :type ci: float
    :return: Cutoff value and its derivative.
    :rtype: float, float
    """
    return 1, 0


@njit
def quadratic_cutoff(r_cut, ri, ci):
    """A quadratic cutoff that goes to zero smoothly at the cutoff boundary.

    :param r_cut: Cutoff value (in angstrom).
    :type r_cut: float
    :param ri: Interatomic distance.
    :type ri: float
    :param ci: Cartesian coordinate divided by the distance.
    :type ci: float
    :return: Cutoff value and its derivative.
    :rtype: float, float
    """
    rdiff = r_cut - ri
    fi = rdiff * rdiff
    fdi = 2 * rdiff * ci

    return fi, fdi


@njit
def cubic_cutoff(r_cut, ri, ci):
    """A cubic cutoff that goes to zero smoothly at the cutoff boundary.

    :param r_cut: Cutoff value (in angstrom).
    :type r_cut: float
    :param ri: Interatomic distance.
    :type ri: float
    :param ci: Cartesian coordinate divided by the distance.
    :type ci: float
    :return: Cutoff value and its derivative.
    :rtype: float, float
    """

    rdiff = r_cut - ri
    fi = rdiff * rdiff * rdiff
    fdi = 3 * rdiff * rdiff * ci

    return fi, fdi


@njit
def cosine_cutoff(r_cut, ri, ci, d=1):
    """A cosine cutoff that returns 1 up to
    .. .. math::
    ..    r_{cut} - d,
    and assigns a cosine envelope to values of r between (r_cut - d) and r_cut.

    :param r_cut: Cutoff value (in angstrom).
    :type r_cut: float
    :param ri: Interatomic distance.
    :type ri: float
    :param ci: Cartesian coordinate divided by the distance.
    :type ci: float
    :return: Cutoff value and its derivative.
    :rtype: float, float
    """

    if ri > r_cut - d:
        fi = (1/2) * (cos(pi * (ri - r_cut + d) / d) + 1)
        fdi = (pi/(2 * d)) * sin(pi * (r_cut - ri) / d) * ci
    else:
        fi = 1
        fdi = 0

    return fi, fdi
