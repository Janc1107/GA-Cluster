from matplotlib.pyplot import *
from scipy.cluster.hierarchy import dendrogram, set_link_color_palette
from ase.data.clustering_v3 import AgglomerativeClustering
from ase.data.comparators_opt import BagOfBonds
from ase.io import read,write
from random import choice, sample
from math import ceil
import numpy as np

def pop_filter():
    #data = read('less_candidates.traj@:')
    data = read('candidates.traj@:')
    n_samples = len(data)

# Comparator
    print ('comp')
    comp = BagOfBonds(excluded_types=[1])

# Cluster data
    print ('AC')
    ac = AgglomerativeClustering(comp=comp, linkage='average', unique_data=False)
    feature_matrix, similarity_matrix, linkage_matrix = ac.grow_tree(data=data)

# Cut tree
    test = 1.0
    print ('\n\ninitial test = ', test)

    #print 'ac.cut_tree'
    labels, clusters, branches, centroids, cluster_energies, avg_width = \
        ac.cut_tree(t= test, criterion='inconsistent', cluster_min=2)

#    print 'labels = ', labels
    print ('# of data = ', n_samples)
    print ('# of clusters = ', len(clusters))
    print ('# of branches = ', len(branches))
#    print 'centroids = ', centroids
#    print 'cluster_energies = ', cluster_energies
#    print 'avg_width = ', avg_width

    while len(clusters) > 10:
        print ('\nclusters>10')
        test += 0.1
        print ('test = ', test)

        #print 'ac.cut tree'
        labels, clusters, branches, centroids, cluster_energies, avg_width = \
            ac.cut_tree(t= test, criterion='inconsistent', cluster_min=2)
        print ('clusters:', len(clusters))

    while len(clusters) < 10:
        print ('\nclusters<10')
        test -= 0.1
        print ('test = ', test)
        #print 'ac.cut tree'
        labels, clusters, branches, centroids, cluster_energies, avg_width = \
            ac.cut_tree(t= test, criterion='inconsistent', cluster_min=2)
        print ('clusters:', len(clusters))
        if len(clusters) == 1:
            break

    print ('t-value is:', test)
    n_clusters = len(branches[0])

# Make every cluster an own traj file
    for i in range(len(clusters)):
        starting_population = []
        print ('cluster:', i, 'length = ', len(clusters[i]) )  # Show which cluster is using
        len_cluster = len(clusters[i])

        # If there are more than 20 structures in a cluster, pick 20
        # Else, put all structures in the file
        if clusters[i] != []:
            if len_cluster >= 20:
                starting_population = sample(clusters[i], 20)
            else:
                for j in range(len_cluster):
                        starting_population.append(clusters[i][j])

        print ('length of starting_population = ', len(starting_population))
        print ('branches = ', n_clusters)

        # Output file
        write('new_cluster_{0}.traj'.format(i), starting_population)

#    write('starting_pop_filter.traj',starting_population)
#    return starting_population


"""
## FIGURE ##
    params={'lines.linewidth':1.5,
            'legend.fontsize':8,
            'xtick.labelsize':8,
            'ytick.labelsize':8,
            'axes.labelsize':8,
            'axes.linewidth':0.5}

    rcParams.update(params)
"""
"""
# Get colors for dendrogram branches
    color_palette = ['#a50f15','#084594']
    set_link_color_palette(color_palette)
    n_colors = len(color_palette)
    color_palette = ['grey']+color_palette*int(ceil(n_clusters/float(n_colors)))

    colors = ['grey']*(2*n_samples-1)
    for cluster, label in zip(branches[0],branches[1]):
        twigs = ac.get_leaves(cluster)[1]
        for twig in twigs:
            colors[int(twig)] = color_palette[label]
"""
"""
    starting_population = []
    for i in range(len(clusters)):
        print 'clusters:', clusters[i]
        if clusters[i] != []:
            choise = choice(clusters[i])
            starting_population.append(choise)


    write('starting_pop_filter.traj',starting_population)
    return starting_population
"""
