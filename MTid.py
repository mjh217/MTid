#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
"""
Assess the probability that a seismic moment tensor belongs to each of the known event-
type populations using scaled von Mises-Fisher probability distribution functions and the 
procedure outlined in Hoggard et al. (2024, GJI), which builds on Ford et al. (2020, GJI).

usage: MTid.py [-h] [-v] [MTs.csv] 

input:
  MTs.csv   Comma-separated file containing moment tensor elements ordered according to
            --> Mxx, Mxy, Mxz, Myy, Myz, Mzz
            where elements can be absolute or normalised and in any consistent unit (e.g., N-m, dyne-cm, Ak).

optional arguments:
  -h        Prints a help message.
  -v        Prints the version number.

output: 
  lune_coords.csv  Lune coordinates for the input MTs.
  prob.csv         Comma-separated list of probability of belonging to each population.


version: 2023.12
"""

import numpy as np
import math
import argparse
from numpy.linalg import norm
from scipy.special import iv

#----------------------------------------USER-PARAMETERS---------------------------------------------

# dimensions of data in real space
dim = 6 

def generate_populations():

   population_1 = population(
      event = "explosions",
      prior = 0.9999,
      gamma = np.array([
      [ 0.36740822,  0.28912206,  0.07658822,  0.36036246,  0.28245966,  0.75226945],
      [ 0.41243448,  0.37831385,  0.09172119,  0.57011912, -0.13295873, -0.57935293],
      [ 0.83208951, -0.28548006, -0.11076200, -0.45832301, -0.02211478, -0.05754047],
      [-0.00453744, -0.43562773, -0.77987658,  0.44715794,  0.04086863,  0.01949217],
      [ 0.03888890, -0.70802624,  0.60379566,  0.36350937, -0.00777659,  0.0204391 ],
      [-0.03167187, -0.02674923,  0.02601096, -0.05434901,  0.94885127, -0.30713557]
      ]),
      # optimal SvMF parameters
      kappa = 54.65367804,
      a = np.array([2.0725349917445777, 2.0527225979792427, 1.2762063394892806, 0.47857163522325175, 0.3848575516754583]),
      # optimal vMF parameters (turn on next two to use a vMF instead)
      #kappa = 34.49525335,
      #a = np.array([1, 1, 1, 1, 1]),
      )

   population_2 = population(
      event = "earthquakes",
      prior = 0.0001,
      gamma = np.array([
      [-0.42034001,  0.13441548,  0.80750027,  0.29492754,  0.25714967,  0.00904755],
      [ 0.88547429,  0.22260801,  0.30730861,  0.18300082,  0.16013458, -0.11318128],
      [-0.04782012, -0.01290958, -0.36344511,  0.92298289,  0.00720227,  0.11621081],
      [ 0.10057037,  0.04350852,  0.02754056, -0.10749762,  0.14407665,  0.97720024],
      [-0.06296248, -0.00478161, -0.29012272, -0.10784493,  0.93907639, -0.13544993],
      [-0.15127501,  0.96451835, -0.19100448, -0.06666910, -0.07454231, -0.01833563]
      ]),
      # optimal SvMF parameters
      kappa = 2.60879588,
      a = np.array([1.5258071357847298, 1.0744906950919515, 0.9805314488231607, 0.8035416739833638, 0.7741547209131909]),
      # optimal vMF parameters (turn on next two to use a vMF instead)
      #kappa = 2.54531175,
      #a = np.array([1, 1, 1, 1, 1]),
      )
   
   populations = [population_1, population_2]
   
   # Check the priors sum to 1
   if np.sum([pop.prior for pop in populations]) != 1. :
      print('Error: Priors currently sum to', np.sum([pop.prior for pop in populations]), 'but must sum to 1. Fix inside function generate_populations(). Exiting.')
      exit()

   return populations





def main():

   global n, i

   parser = argparse.ArgumentParser(
      description=('Assess the probability that a seismic moment tensor belongs to each of the known event-type '
      'populations using scaled von Mises-Fisher probability distribution functions and the procedure outlined in '
      'Hoggard et al. (2024, GJI), which builds on Ford et al. (2020, GJI).'),
      formatter_class=argparse.RawTextHelpFormatter)
   parser.add_argument('MTs', help='Comma-separated file containing six moment tensor elements ordered according to  \n'
      '--> Mxx, Mxy, Mxz, Myy, Myz, Mzz  \n'
      'where elements can be absolute or normalised and in any consistent unit (e.g., N-m, dyne-cm, Ak).')
   parser.add_argument('-v','--version', action='version', version='2023.12')
   args = parser.parse_args()
   print("Reading input MTs from file:",args.MTs)

   # Get the statistics for each population
   populations = generate_populations()

   # Read input MTs into an array
   MTs = np.loadtxt(fname=args.MTs, delimiter=',', dtype='float')
   n = np.shape(MTs)[0]
   print("Found",n,"MT solutions. Calculating probability of belonging to each one of", len(populations),"event-type populations...")

   # Normalise and reproject to six-degree unit vector
   x = np.zeros((n, dim))
   for j in range(0, n):
	   x[j,:] = reproject_MT(MTs[j,:])

   # Loop through the MTs and find population probabilities
   f_test = np.zeros((n,len(populations)))
   P = np.zeros((n,len(populations)))

   # Evaluate the SvMF for each population 
   f = np.array([[pop.evaluate_SvMF_distribution(x[j, :]) for pop in populations] for j in range(n)])

   # Evaluate the probability of belonging to each population
   priors = np.array([[pop.prior for pop in populations]]*n)
   P = np.einsum("i, ij->ij", 1/np.sum(f*priors, axis=1), f*priors)

   # Raise a warning if the MT most likely belongs to the first group
   first_population = populations[0]
   for j in range(0, n):
      if any (value < P[j,0] for value in P[j,1:]):
         print(f"Event {j+1} likely belongs to {first_population.event}.")

   # Write PDF values and probabilities to outfile
   with open("prob.csv", mode='w') as outfile:
      outfile.write("\n".join([",".join([f"{j}" for j in i]) for i in P]))
      
   # Project MTs onto the fundamental lune of Tape & Tape (2012) for plotting
   lune_coords = np.zeros((n, 2))
   for j in range(0, n):
	   lune_coords[j,:] = lune_MT(MTs[j,:])
   with open("lune_coords.csv", mode='w') as outfile:
      outfile.write("\n".join([",".join([f"{j}" for j in i]) for i in lune_coords]))



class population(object):
   def __init__(self, event, prior, gamma, kappa, a):
      self.event = event
      self.prior = prior
      self.gamma = gamma
      self.kappa = kappa
      self.a = a
   
   def evaluate_SvMF_distribution(self, x):
      """
      Evaluate SvMF function for this population at this location (Equation 9).
      """
      global gamma, kappa, a
      gamma = self.gamma
      kappa = self.kappa
      a = self.a
      return evaluate_SvMF_distribution(x)
      


def reproject_MT(MTs):
   """
   Normalise MT and map to its six-degree unit vector (Equation 1).
   """
   x = np.array([MTs[0], MTs[3], MTs[5], MTs[4]*math.sqrt(2), MTs[2]*math.sqrt(2), MTs[1]*math.sqrt(2)])
   return x/norm(x)



def lune_MT(MT):
   """
   Project MT onto fundamental lune of Tape & Tape (2012)
   """
   # construct symmetric MT matrix M and normalise
   M = [[MT[0], MT[1], MT[2]], [MT[1], MT[3], MT[4]], [MT[2], MT[4], MT[5]]]
   Mnorm = M/np.linalg.norm(M)
   # extract eigenvalues and convert to lune coordinates
   eigval, eigvec = np.linalg.eig(Mnorm)
   eigval_sort = np.sort(eigval)
   gamma = (180/math.pi)*math.atan( (- eigval_sort[2] + 2*eigval_sort[1] - eigval_sort[0]) / (math.sqrt(3)*(eigval_sort[2]-eigval_sort[0])))
   delta = (180/math.pi)*(0.5*math.pi - math.acos( (eigval_sort[2] + eigval_sort[1] + eigval_sort[0]) / math.sqrt(3) ))
   return gamma, delta



def evaluate_SvMF_distribution(x):
   """
   Evaluate SvMF function at this location (Equation 9).
   """
   # Check that the dimensions of data and population parameters are consistent
   if (x.size != dim):
      raise Exception("Data do not possess correct number of dimensions (",dim,")... Exiting.")
   if (np.shape(gamma)[0] != dim or np.shape(gamma)[1] != dim):
      raise Exception("Gamma parameter does not have correct",dim,"x",dim,"entries... Exiting.")
   if (a.size != dim-1):
      raise Exception("Shape parameter [a] does not possess correct number of entries (",dim-1,")... Exiting.")

   # Evaluate the prefactor including calling the modified Bessel function of the first kind
   prefactor = (kappa**((dim/2.)-1.)) / (((2.*np.pi)**(dim/2.)) * iv(((dim/2.) - 1.), kappa))

   # Add a leading 1 to the shape parameter and evaluate the shape factor summation term
   a_extended = np.hstack(([1], a))
   c = np.sum((np.dot(x,gamma)/a_extended)**2)
   
   # Evaluate the exponential term
   exp_factor = np.exp(kappa * np.dot(x, gamma[:,0]) / c**0.5 )

   return prefactor * c**(-0.5*(dim-1.)) * exp_factor



if __name__ == '__main__':
   main()
