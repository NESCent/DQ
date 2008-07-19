#!/usr/bin/env python

# Using CapWords for function, method, and class names
# Using underscored_names for variable names
# Using firstWordLower for module, package names.
# Using ALL_CAPS_WITH_UNDERSCORES for file handles

#standard imports 
import sys
import os
import re
#from itertools import *

#python extensions
#from numpy import * 

#Ganesh's modules
from ldComputations import *

def RandomIndices(n, k):
    """ return k unique indices at random from range(0, n): i.e., [0, n-1]
    (both included) """

    r_indices = []
    while len(r_indices) < k:
        candidate = random.randint(0, n)
        if candidate in r_indices:
            continue
        else:
            r_indices.append(candidate)
    r_indices.sort()        
    print 'INDICES = ', r_indices
    return r_indices

def HashIndices(index_list, tolerance):

    index_list.sort()
    k = len(index_list)
    distances = [index_list[i] - index_list[i-1] for i in range(1, k)]

    mod_what = tolerance + 1
    remainders = [distances[i]%mod_what for i in range(0,len(distances))] 

    binned_distances = [distances[i] - remainders[i] for i in range(0, len(distances))]
    #print 'BINS = ', binned_distances
    
    if (k == 3):
        return sorted(binned_distances)
    else:
        return binned_distances

def Factorial(n):
    fact = 1
    for i in range(1, n+1):
        fact = fact*i
    return fact    


def NchooseK(n, k):
    return Factorial(n)/(Factorial(n-k) * Factorial(k))


def AppendExtend(a, b):
    w = None
    if isinstance(a, list) and isinstance(b, list):
        w = a + b
        #print 'In Append 11', a, b, w

    if isinstance(a, list) and not(isinstance(b, list)):
        w = a + [b]
        #print 'In Append 10', a, b, w

    if not isinstance(a, list) and isinstance(b, list):
        w = [a] + b
        #print 'In Append 01', a, b, w

    if not isinstance(a, list) and not isinstance(b, list):
        w = [a, b]
        #print 'In Append 00', a, b, w
    return w
        
def CartesianProduct(list1, list2):
    """ the "cartesian product" of the two lists. But the lists are 
    not treated as sets: for example, [1, 1] x [3, 2] = [[1, 3], [1, 2],
    [1, 3], [1, 2]]; further, [a, b] x [c] = [a.extend(c), b.extend(c) if
    a, b and c are lists"""

    product = [AppendExtend(i, j) for i in list1 for j in list2] #as simple as that!!!
    #print 'HERE ', list1, list2
    return product


def ExponentialDistances(n):
    """return [n/2, (n/2)/2, ((n/2)/2)/2..., 1]"""

    #k = n/2
    k = (n-1)/2
    distances = []
    while k > 0:
        distances.append(k)
        k = k/2
    #print 'DISTANCES = ', distances    
    return distances


def MultilocusExponentialDistances(n, level):
    """level here is the number loci involved in the diseq
    computatation. n is the total number of loci. 
    
    Each set of distances is
    a level-1 dimensional vector. Along each dimension, the permitted
    values are (as in ExponentialDistances) (n-1)/2, (n-1)/2/2, ..., 1
    Thus there are about log(n) possible values in each dimension. 
    Further, for a vector (d_1, d_2, d_3, ..., d_{level-1}),
    d_1+d_2+...+d_{level-1} <= n-1

    The idea is that given the above vector of distances, we look at all
    partitions (i.e., set of loci) l_1, l_2, ..., l_{level} such that
    l_2-l_1 = d_1, l_3-l_2 = d_2 etc.
    """
    

    along_one_dimension = ExponentialDistances(n)
    along_all_dimensions = along_one_dimension
    for i in range(2, level):
        #print 'ALONGS = ', along_all_dimensions, along_one_dimension
        along_all_dimensions = CartesianProduct(along_all_dimensions, along_one_dimension)
    
    #print 'ALONGS = ', along_all_dimensions, along_one_dimension
    #amazing list comprehension!!!
    filtered = [item for item in along_all_dimensions if sum(item) < n]
    return filtered            

def PartitionsGivenDistanceVector(n, distance_vector):
    """given distance_vector = [d_1, d_2, ..., d_k], and the total #loci =
    n; return one by one partitions (l1, l2,..., l_{k+1}) such that l_i in range(0, n) and
    l_{i+1} - l_i = d_i; this is an iterator; so after all the partitions
    have been yielded, keep returning None"""

    # ugly special case dealing with when distance_vector may just be a
    # scalar distance
    if not isinstance(distance_vector, list):
        distance_vector = [distance_vector]

    end_to_end_distance = sum(distance_vector)
    cumsum = [sum(distance_vector[0:k+1]) for k in range(0, len(distance_vector))]
        
    for i in range(0, n-end_to_end_distance):
        yield [i] + [i + d for d in cumsum]

    while True:
        yield None


def AllPartitionsGivenLevel(n, level):
    """level is the #loci we want in each partition. If level = 3, we want
    three-locus partitions. But we are *not* generating all the \theta(n^3)
    partitions."""

    all_distance_vectors = MultilocusExponentialDistances(n, level)
    for d_vector in all_distance_vectors:
        all_partitions = PartitionsGivenDistanceVector(n, d_vector)
        partition = all_partitions.next()
        while not partition == None:
            yield partition
            partition = all_partitions.next()
    while True:
        yield None

    

def AllPairsOfLoci(n):
    """return [(i, j) such that i, j in range(0, n), and abs(i-j) in [n/2,
       (n/2)/2, ((n/2)/2)/2, ..., 1]]. For n = 10, the lisy
       but this is function is an iterator; so it yields the pairs in the
       above set one by one until all pairs have been yielded, and then keeps
       yielding None"""
    
    exponential_distances = ExponentialDistances(n)
    #print exponential_distances
    for distance in exponential_distances:
        for i in range(0, n-distance):
            yield [i, i+distance]
    while True:    
        yield None    

def PairsOfLociGivenDistance(n, distance):
    """return [(i, j) such that i, j in range(0, n), and abs(i-j) = distance]
       but this is function is an iterator; so it yields the pairs in the
       above set one by one until all pairs have been yielded, and then keeps
       yielding None"""

    for i in range(0, n-distance):
        yield [i, i+distance]
    while True:    
        yield None    

def OnePartition(haplotype_matrix, partition):
    (n_rows, n_col) = haplotype_matrix.shape
    partitioned_matrix = haplotype_matrix[0:n_rows, partition]
    z = PartitionDisEq(partitioned_matrix, partition)
    return z

def OneFixedDist(haplotype_matrix, distance):
    (n_rows, n_col) = haplotype_matrix.shape
    print n_rows, n_col

    #n_col = 10
    #pairs = PairsOfLociGivenDistance(n_col, distance)
    partitions = PartitionsGivenDistanceVector(n_col, distance)
    lds = [] 

    #partition = pairs.next()
    partition = partitions.next()
    while not(partition == None):
        print partition
        z = OnePartition(haplotype_matrix, partition)
        lds.append(z)
        partition = partitions.next()
    return lds        


def AllFixedDistances(haplotype_matrix, level, tolerance):

    (n_rows, n_col) = haplotype_matrix.shape
    print n_rows, n_col

    #n_col = 10
    #all_pairs = AllPairsOfLoci(n_col)
    all_partitions = AllPartitionsGivenLevel(n_col, level)
    fixed_distances_lds = [] 
    distance_lists = []

    partition = all_partitions.next()
    while not(partition == None):
        print partition
        z = OnePartition(haplotype_matrix, partition)
        distance_bin = HashIndices(partition, tolerance)
        #print distance_bin, partition
        if not(distance_bin in distance_lists):
            distance_lists.append(distance_bin)
            fixed_distances_lds.append([z])
        else:
            hash_value = distance_lists.index(distance_bin)
            fixed_distances_lds[hash_value].append(z)
        partition = all_partitions.next()

    return [distance_lists, fixed_distances_lds]                 
