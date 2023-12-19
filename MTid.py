#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
import argparse
from collections import namedtuple

import numpy as np
from scipy.special import iv

from populations import a, dim, events, gamma, kappa, prior


def parser():
    parser = argparse.ArgumentParser(
        description='Moment Tensor identification \n \n'
                    'Assess the probability that a seismic moment tensor belongs to each of the known event-type '
                    'populations using scaled von Mises-Fisher probability distribution functions and the '
                    'procedure outlined in Hoggard et al. (2024, GJI), which builds on Ford et al. (2020, GJI).',
        epilog='output: \n'
               '  lune_coords.csv  Lune coordinates for the input MTs. \n'
               '  prob.csv         Comma-separated list of probability of belonging to each population.',
        formatter_class=argparse.RawTextHelpFormatter, )
    parser.add_argument("MTs",
        help='Comma-separated file containing six moment tensor elements ordered according to\n'
             '--> Mxx, Mxy, Mxz, Myy, Myz, Mzz\n'
             'where elements can be absolute or normalised and in any consistent unit (e.g., N-m, dyne-cm, Ak).',
        metavar="csv_file", )
    parser.add_argument("-v", "--version", action="version", version="%(prog)s 2023.12")

    return parser.parse_args()


def get_PDF_parameters(events, prior, gamma, kappa, a):
    """
    Read the parameters associated with each population.
    """
    # Assign the population parameters
    Population = namedtuple("Population", ["event", "prior", "gamma", "kappa", "a"])
    populations = []
    input_vars = [prior, gamma, kappa, a]
    for event in events:
        populations.append( Population(event, *(input_var[event] for input_var in input_vars)) )
    # Perform several checks that the population parameters have the correct dimensions
    assert sum(pop.prior for pop in populations) == 1, "Priors currently sum to {sum(pop.prior for pop in populations)} but must sum to 1."
    for pop, array in gamma.items():
        assert all(shape == dim for shape in array.shape), (f"Provided shape of gamma array for {pop} is {array.shape} but must have dimensions ({dim}, {dim}).")
    for pop, array in a.items():
        assert all(shape == dim - 1 for shape in array.shape), (f"Provided length of shape vector (a) for {pop} is {array.shape[0]} but must have {dim - 1} entries.")

    return populations


def evaluate_SvMF_distribution(MTs, pop):
    """
    Evaluate SvMF function at this location (Equation 9).
    """
    # Evaluate the prefactor including calling the modified Bessel function of the first kind
    prefactor = ( pop.kappa ** ((dim / 2) - 1) / (2 * np.pi) ** (dim / 2) / iv((dim / 2) - 1, pop.kappa) )
    # Add a leading 1 to the shape parameter and evaluate the shape factor summation term
    c = np.sum((np.dot(MTs, pop.gamma) / np.concatenate(([1], pop.a))) ** 2, axis=1)
    # Evaluate the exponential term
    exp_factor = np.exp(pop.kappa * np.dot(MTs, pop.gamma[:, 0]) / np.sqrt(c))

    return prefactor * c ** (-0.5 * (dim - 1)) * exp_factor


def lune_MT(MT):
    """
    Project MT onto fundamental lune of Tape & Tape (2012).
    """
    # Construct symmetric MT matrix M and normalise
    M = [[MT[0], MT[1], MT[2]], [MT[1], MT[3], MT[4]], [MT[2], MT[4], MT[5]]]
    M_norm = M / np.linalg.norm(M)
    # Extract eigenvalues and convert to lune coordinates
    eigvals = np.linalg.eigvalsh(M_norm)
    gamma = np.rad2deg(np.arctan((-eigvals[2] + 2 * eigvals[1] - eigvals[0]) / (np.sqrt(3) * (eigvals[2] - eigvals[0]))))
    delta = np.rad2deg( np.pi / 2 - np.arccos((eigvals[2] + eigvals[1] + eigvals[0]) / np.sqrt(3)) )

    return gamma, delta


def main():

    # Read in command line arguments and extract PDF parameters
    args = parser()
    populations = get_PDF_parameters(events, prior, gamma, kappa, a)

    # Read input MTs into an array
    print("Reading input MTs from file:", args.MTs)
    MTs = np.loadtxt(fname=args.MTs, delimiter=",", dtype="float")
    assert MTs.shape[1] == dim, "Number of elements does not match {dim} for at least one moment tensor."
    n = MTs.shape[0]
    print(f"Found {n} MT solutions. Calculating probability of belonging to each one of {len(populations)} event-type populations...")

    # Normalise and reproject to six-degree unit vector (Equation 1)
    MTs_proj = MTs[:, [0, 3, 5, 4, 2, 1]]
    MTs_proj[:, -3:] *= np.sqrt(2)
    MTs_proj = MTs_proj / np.linalg.norm(MTs_proj, axis=1).reshape(-1, 1)

    # Evaluate the SvMF for each population (Equation 9)
    f = np.column_stack([evaluate_SvMF_distribution(MTs_proj, pop) for pop in populations])

    # Evaluate the probability of belonging to each population (Equation 17) and write out
    priors = [pop.prior for pop in populations]
    P = np.einsum("i, ij->ij", 1 / np.sum(f * priors, axis=1), f * priors)
    np.savetxt("prob.csv", P, delimiter=",")

    # Raise a warning if the MT most likely belongs to the first group
    for event_number in np.where(np.argsort(P, axis=1)[:, 0] == len(populations) - 1)[0]:
        print(f"Event {event_number + 1} likely belongs to {populations[0].event}.")

    # Project MTs onto the fundamental lune of Tape & Tape (2012) for plotting purposes
    lune_coords = np.zeros((n, 2))
    for i, MTs_row in enumerate(MTs):
        lune_coords[i] = lune_MT(MTs_row)
    np.savetxt("lune_coords.csv", lune_coords, delimiter=",", fmt="%7.3f")


if __name__ == "__main__":
    main()
