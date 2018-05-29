# Main.py
# This is the main function to run in the project
#from ase.ga.data import PrepareDB, DataConnection
from molemaker import Starting_population
#from random import random
from ase.io import write
#from ase.optimize import BFGS
#from ase.calculators.dftb import Dftb
from ase import Atoms
import numpy as np
#from ase.ga.population import Population
#from ase.ga.standard_comparators import InteratomicDistanceComparator
#from ase.ga.cutandsplicepairing import CutAndSplicePairing
#from ase.ga.utilities import closest_distances_generator
#from ase.ga.utilities import get_all_atom_types
#from ase.ga.offspring_creator import OperationSelector
#from ase.ga.standardmutations import MirrorMutation
#from ase.ga.standardmutations import RattleMutation
#from ase.ga.standardmutations import PermutationMutation
from ase.constraints import FixAtoms
#from starting_pop_maker import pop_maker
from starting_pop_filter import pop_filter
from GA import GA
from time import time
import platform
from os.path import exists
# New
import multiprocessing as mp


#print platform.node()
"""
# create the surface
slab = Atoms('', positions = np.zeros((0,3)), cell = [30., 30., 30.])
slab.set_constraint(FixAtoms(mask=len(slab) * [True]))

# Define the volume in which the adsorbed cluster is optimized
# the volume is defined by a corner position (p0)
# and three spanning vectors (v1, v2, v3)
pos = slab.get_positions()
cell = slab.get_cell()
p0 = np.array([0., 0., 12]) #max(pos[:, 2]) + 2.])
v1 = cell[0, :] * 0.8
v2 = cell[1, :] * 0.8
v3 = cell[2, :]
v3[2] = 3.
"""

# Define the composition of the atoms to optimize
atom_numbers = 9 * [6] + 7 * [1] + 1 * [7]
#atcount = 5000
atcount = 500
moleculegroup = []

# Make initial candidates
start = time()
ga = Starting_population(atom_numbers, moleculegroup)
ga.randomizer()
while atcount > len(moleculegroup):
    ga.randomizer()
end = time()
print ('time to generate starting population:', end - start)

write('candidates.traj', moleculegroup)

# Cluster
# Search in a group
start = time()
#starting_filter = pop_filter()
pop_filter()
end = time()
print ('time to cluster & create starting population:', end - start)

# Do genetic algorithm on every cluster
# Do open MP
for cluster in range(20):
        if isfile('cluster_{0}.traj'.format(cluster)) == 0 :    # If the doesn't file exist
                #print 'cluster_{0}.traj'.format(i)
                #GA(cluster = i, atom_numbers = atom_numbers)
                print 'There is no more cluster'
                break
        else:                                           # If the file exist
                #GA(cluster=cluster, atom_numbers=atom_numbers)
                pool = mp.Pool()
                args = [cluster, atom_numbers]
                result = pool.map(target = GA, args)