# This file is used to do GA

from ase.ga.data import PrepareDB, DataConnection
#from molemaker import Starting_population
from random import random
from ase.io import write
from ase.optimize import BFGS
from ase.calculators.dftb import Dftb
from ase import Atoms
import numpy as np
from ase.ga.population import Population
from ase.ga.standard_comparators import InteratomicDistanceComparator
from ase.ga.cutandsplicepairing import CutAndSplicePairing
from ase.ga.utilities import closest_distances_generator
from ase.ga.utilities import get_all_atom_types
from ase.ga.offspring_creator import OperationSelector
from ase.ga.standardmutations import MirrorMutation
from ase.ga.standardmutations import RattleMutation
from ase.ga.standardmutations import PermutationMutation
from ase.constraints import FixAtoms
#from starting_pop_maker import pop_maker
from starting_pop_filter import pop_filter
from time import time
#import platform


def GA(cluster, atom_numbers):
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

	# Make a db file
	db_file = 'gadb_{0}.db'.format(cluster)

    # Create the database to store information in
    starting_pop = read('cluster_{0}.traj@:'.format(cluster))
    print 'This is cluster', cluster
    if len(starting_pop) <  2:      # If there has no enough structures to do GA, skip that cluster
        print "This cluster doesn't has enough structures to do GA."
        #continue
        #break
        return


	# Modify from here
	# create the database to store information in
	starting_pop = read('cluster_{0}.traj@'.format(cluster))

	d = PrepareDB(db_file_name=db_file,
					simulation_cell = slab,
					stoichiometry=atom_numbers)
	

	for i in starting_pop:
		d.add_unrelaxed_candidate(i)

	population_size = len(starting_pop)
	mutation_probability = 0.0
	n_to_test = 5000  # Number to test new candidates

	da = DataConnection('gadb_{0}.db'.format(cluster))
	atom_numbers_to_optimize = da.get_atom_numbers_to_optimize()
	n_to_optimize = len(atom_numbers_to_optimize)
	slab = da.get_slab()
	all_atom_types = get_all_atom_types(slab, atom_numbers_to_optimize)
	blmin = closest_distances_generator(all_atom_types,
										ratio_of_covalent_radii=0.7)

	comp = InteratomicDistanceComparator(n_top=n_to_optimize,
										pair_cor_cum_diff = 0.015,
										pair_cor_max=0.7,
										dE=0.02,
										mic=False)

	pairing = CutAndSplicePairing(slab, n_to_optimize, blmin)
	mutations = OperationSelector([1., 1., 1.],'/'
								[MirrorMutation(blmin, n_to_optimize),
								RattleMutation(blmin, n_to_optimize),
								PermutationMutation(n_to_optimize)])

	while da.get_number_of_unrelaxed_candidates() > 0:
		start = time()
		a = da.get_an_unrelaxed_candidate()
		a.set_calculator(Dftb(label='C9H7N',
							Hamiltonian_MaxAngularMomentum_='',
							Hamiltonian_MaxAngularMomentum_C='"p"',
							Hamiltonian_MaxAngularMomentum_H='"s"',
							Hamiltonian_MaxAngularMomentum_N='"p"'
							))
		print('Relaxing starting candidate {0}'.format(a.info['confid']))
		dyn = BFGS(a, trajectory=None, logfile=None)
		dyn.run(fmax=0.05, steps=100)
		a.info['key_value_pairs']['raw_score'] = -a.get_potential_energy()
		da.add_relaxed_step(a)
		end = time()
		print 'time to optimize:', end - start

	population = Population(data_connection=da,
						population_size=population_size,
						comparator=comp)

	# test n_to_test new candidates
	for i in range(n_to_test):
		print('Now starting configuration number {0}'.format(i))
		a1, a2 = population.get_two_candidates()
		a3, desc = pairing.get_new_individual([a1, a2])
		if a3 is None:
			continue
		da.add_unrelaxed_candidate(a3, description=desc)

		# Check if we want to do a mutation
		if random() < mutation_probability:
			a3_mut, desc = mutations.get_new_individual([a3])
			if a3_mut is not None:
				da.add_unrelaxed_step(a3_mut, desc)
				a3 = a3_mut

		# Relax the new candidate
		start = time()
		a3.set_calculator(Dftb(label='C9H7N',
							Hamiltonian_MaxAngularMomentum_='',
							Hamiltonian_MaxAngularMomentum_C='"p"',
							Hamiltonian_MaxAngularMomentum_H='"s"',
							Hamiltonian_MaxAngularMomentum_N='"p"'
							))

		dyn = BFGS(a3, trajectory=None, logfile=None)
		dyn.run(fmax=0.05, steps=100)
		a3.info['key_value_pairs']['raw_score'] = -a3.get_potential_energy()
 		da.add_relaxed_step(a3)
		population.update()
		end = time()
		print 'time to optimize:', end - start
		if a3.get_potential_energy() < -557.783851467:
			break

	write('all_candidates_in_cluster{0}.traj'.format(cluster), da.get_all_relaxed_candidates())

