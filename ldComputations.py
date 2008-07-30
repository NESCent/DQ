#!/usr/bin/env python

# Using CapWords for function, method, and class names
# Using underscored_names for variable names
# Using firstWordLower for module, package names.
# Using ALL_CAPS_WITH_UNDERSCORES for file handles

#standard imports 
import sys
import os
import re

#python extensions
from numpy import * 

import levelDistancePartition

############################
############################

def PerHaplotypeScoresVector(haplotype, number_of_ones, number_of_haplotypes):
    """
    Compute and return the vector of disequilibrium scores for the given haplotype.

    If the haplotype is of length N, the vector is of length 2**N - 1.

    Input parameters
    ----------------

    haplotype                       is a 1-dimensional 0-1 numpy array,
                                    
    number_of_ones                  is a 1-dimensional numpy array. 
                                    number_of_ones[j] is
                                    the number of type 1 alleles in the
                                    sample at site j (counting from 0)

    number_of_haplotypes            number of haplotypes in the sample from
                                    which the input came.

    Return value
    ------------

    A vector row_scores[0:2**N-1] (where N is the length of the input haplotype), such that 
    row_scores[j] is the multilocus disequilibrium computed for the partition of haplotype 
    represented by the binary representation of j+1. 
    For example, let N = 4,
    and j = 5. Then, row_scores[j] is the multilocus disequilibrium of the
    partition of haplotype specified by [0, 1, 1, 0] i.e., of the partition
    [haplotype[1], haplotype[2]]. So in this example, row_scores[5] is the
    individual haplotype version of the traditional 2-locus linkage
    disequilibrium between sites 1 & 2.

    To see how a multilocus disequilibrium score is computed for a given
    partition of a single haplotype, see function HaplotypePartitionScore
    """

    n_cols = len(haplotype)
    if n_cols < 2:
        print 'Cannot calculate disequilibria for a single column ', partition
        return -1

    row_scores = []
    all_partitions = levelDistancePartition.AllPartitions(n_cols)
    partition = all_partitions.next()
    while not partition == None:
        #print partition
        score = HaplotypePartitionScore(haplotype, number_of_ones, number_of_haplotypes, partition)
        row_scores.append(score)
        partition = all_partitions.next()
    return row_scores    


def HaplotypePartitionScore(haplotype, number_of_ones, number_of_haplotypes, partition):
    """
    Compute the multi-locus disequilibrium score for a specified
    haplotype partition.
    

    Input parameters
    ----------------

    haplotype                       is a 1-dimensional 0-1 numpy array,
                                    
    number_of_ones                  is a 1-dimensional numpy array. 
                                    number_of_ones[j] is
                                    the number of type 1 alleles in the
                                    sample at site j (counting from 0)

    partition                       is a 0-1, 1-dimensional, non-numpy, array 
                                    representing a selection of sites or
                                    columns.

    number_of_haplotypes            number of haplotypes in the sample from
                                    which the input came.

    Return value & Details
    ----------------------

    returns the multi-locus disequilibrium calculated for the specified haplotype partition. 
    The detailed calculations are below, where the return value is denoted row_diseq.
                                                   =================================

    Let:

    1. haplotype = h[0:N] (i.e., an array of length N, remembering that 0:n = {0, 1, 2, ..., N-1}
    2. partition = p[0:N] 
    3. number_of_haplotypes = Q

    Then, HaplotypePartitionScore's calculations are:

    partitioned_haplotype[]                 [h[j] such that partition[j] == 1] (just those
                                            sites j for which partition[j] == 1)

    partitioned_number_of_ones[]            [number_of_ones[j] for j where partition[j] == 1]

    partitioned_number_of_zeros[]           Q-partitioned_number_of_ones
                                            (subtract each element of partitioned_number_of_ones
                                            from Q)

    product_of_fractions_ones               partitioned_number_of_ones.prod()/(Q**N)
                                            [= Product_{j from 0 through N-1} number_of_ones[j]*partition[j]/Q]

    product_of_fractions_zeros              partitioned_number_of_zeros.prod()/(Q**N)
                                            [= Product_{j from 0 through N-1} number_of_zeros[j]*partition[j]/Q]

    indicator_all_ones                      = 1, if partitioned_haplotype is made of all ones.
                                            = 0, otherwise

    indicator_all_zeros                     = 1, if partitioned_haplotype is made of all zeros.
                                            = 0, otherwise

    row_diseq                               row_diseq is calculated as:                                            

    indicator_all_zeros + indicator_all_ones - product_of_fractions_zeros - product_of_fractions_ones
    -------------------------------------------------------------------------------------------------
                    1 - product_of_fractions_ones - product_of_fractions_zeros
                                            
    """

    n_sites = len(haplotype)

    # get a list of sites represented in partition.
    sites_in_partition = filter(lambda j: partition[j] == 1, range(0, n_sites))
            
    # partitioned_haplotype contains just those elements of haplotype whose
    # indices are in sites_in_partition. And then array() typecasts it as a
    # numpy array.
    partitioned_haplotype = array(haplotype[sites_in_partition])

    partitioned_number_of_ones = array(number_of_ones[sites_in_partition]) 
    partitioned_number_of_zeros = number_of_haplotypes - partitioned_number_of_ones

    #again, these are whole array operations
    partitioned_fraction_of_ones = partitioned_number_of_ones/float(number_of_haplotypes)
    partitioned_fraction_of_zeros = partitioned_number_of_zeros/float(number_of_haplotypes)

    product_of_fractions_ones = partitioned_number_of_ones.prod()
    product_of_fractions_zeros = partitioned_number_of_zeros.prod()
    
    #product_of_fractions_ones = float(partitioned_number_of_ones.prod())/(number_of_haplotypes**n_sites)
    #product_of_fractions_zeros = float(partitioned_number_of_zeros.prod())/(number_of_haplotypes**n_sites)

    indicator_all_ones = (0 not in partitioned_haplotype)
    indicator_all_zeros = (1 not in partitioned_haplotype)

    nr = indicator_all_zeros + indicator_all_ones - product_of_fractions_zeros - product_of_fractions_ones
    dr = 1 - product_of_fractions_ones - product_of_fractions_zeros

    row_diseq = nr/dr
    return row_diseq


def BitArray(d, n_bits):
    """Puts the n-bit binary representation of positive integer d in a
    numpy array and returns the array
    """

    #x is a numpy array
    x = zeros(n_bits)
    for j in range(0, n_bits):
        x[n_bits-j-1] = fmod(d/pow(2, j), 2)
    return x

############################
############################

def FilterMissingData(m):
    """
    Returns a matrix without the rows in the input matrix m that contain
    -1s

    Input parameter
    ---------------

    m                       a primarily 0-1 *numpy* matrix, with -1 
                                            =======
                            missing data.

    Return value
    ------------

    A 0-1 numpy matrix without the rows in m that contain -1s.

    """

    # filter rows with -1 (missing data) (note: this doesn't modify m)
    matrix_without_missing_data = filter(lambda row: -1 not in row, m)
    # the statement above returns a normal list of numpy arrays. The following
    # statement typecasts the normal list into a numpy array of numpy
    # arrays, like partition_matrix is. Note that we are not modifying
    # partition_matrix here at all.
    matrix_without_missing_data = array(matrix_without_missing_data)

    return(matrix_without_missing_data)



def AlleleOneSiteFrequencies(matrix):
    """
    Computes the allele-one frequency for each site (column) in matrix 
    
    Input parameters
    ----------------

    matrix                  a purely 0-1 *numpy* matrix
                                         =======

   
    Return value
    ------------

    allele_one_site_frequencies[]       allele_one_site_frequencies[j] = number of ones in the
                                        j-th column (counting from 0) of matrix.

    """


    # each row in matrix.transpose() is a column in
    # matrix. note: calling the transpose
    # method of the array object does not modify the object itself.
    allele_one_site_frequencies = array([row.sum() for row in matrix.transpose()])
    return(allele_one_site_frequencies)


def PartitionDisEq(partition_matrix, partition):
    """
    Computes a multi-locus disequilibrium measure for a given matrix
    (partition_matrix).

    The diseq. is computed is E, as defined at the bottom of page 2 in
    Dpca.pdf

    Input parameters
    ----------------

    partition_matrix        a primarily 0-1 *numpy* matrix, with -1 also being used to
                                            =======
                            denoted missing data. ***Note: partition_matrix
                            is treated as a haplotype matrix would be.
                            Specifically, PartitionDisEq does not treat
                            successive rows as representing the
                            maternally-derived and paternal-derived
                            chromosomes of a diploid zygote.***
    
    partition               well, partition_matrix is usually a sub-matrix
                            of some matrix in the calling function that has
                            the same number of rows as the original matrix,
                            but whose columns are a subset of the columns
                            of the original matrix. partition is a list the
                            columns of the original matrix that
                            partition_matrix represents. For example, if
                            partition = [1, 4, 9, 13], it means
                            partition_matrix contains columns 1, 4, 9, 13
                            (in that order) of the original matrix in the 
                            calling function.

    Return value
    ------------

    Let H[M, N] be the input M x N matrix. PartitionDisEq disregards any
    rows that contain a -1 (missing data). So let H[Q, N] be the 0-1 matrix
    with no missing data, where Q is the number of rows with no missing
    data.

    Then, PartitionDisEq's calculations are:

    number_of_ones[i]           number of 1 alleles in column i, i in {0, 1, 2, ..., N-1}
    fraction_ones[i]            number_of_ones[i]/Q
    number_of_zeros[i]          Q-number_of_ones[i]
    fraction_zeros[i]           number_of_zeros[i]/Q = 1-fraction_ones[i]

    product_of_fractions_ones   fraction_ones.prod()  [In practice this may
                                                       be better calculated as
                                                       number_of_ones.prod()/Q**N,
                                                       to avoid floating
                                                       point complications.]
    
    product_of_fractions_zeros  fraction_zeros.prod() [similar remark applies]
    
    indicator_all_one_row[]     indicator_all_one_row[i] = 1 if row_i is
                                all ones, and zero otherwise.
    indicator_all_zero_row[]    indicator_all_zero_row[i] = 1 if row_i is
                                all zeros, and zero otherwise.


    row_diseq[]                 row_diseq[i] is calculated as:

    indicator_all_zero_row[i] + indicator_all_one_row[i] - product_of_fractions_zeros - product_of_fractions_ones
    -------------------------------------------------------------------------------------------------------------
                    1 - product_of_fractions_ones - product_of_fractions_zeros


    matrix_diseq               row_diseq.sum()/Q

    The return value is matrix_diseq

    Details
    -------

    The earlier version of PartitionDisEq calculated the same final return
    value, but ordered the computations differently. But the new way of
    organizing (specifically, calculating matrix_diseq as
    row_diseq.sum()/Q) has the advantage that it cleanly expresses
    matrix_diseq as the expected value of Marcy's individual haplotype
    scores (see Section 2, Dpca.pdf) 
    """

    (n_rows, n_cols) = partition_matrix.shape

    if n_cols < 2:
        print 'Cannot calculate disequilibrium for the single column ', partition
        return -1

    matrix_without_missing_data = FilterMissingData(partition_matrix)
    (n_rows, n_cols) = matrix_without_missing_data.shape
    number_of_ones = AlleleOneSiteFrequencies(matrix_without_missing_data)
    
    # BEGIN block comments
    #
    # number_of_zeros = n_rows - number_of_ones
    # product_of_fractions_ones  = float(number_of_ones.prod())/(n_rows**n_cols)
    # product_of_fractions_zeros = float(number_of_zeros.prod())/(n_rows**n_cols)

    # indicator_all_one_row = [0 not in row for row in matrix_without_missing_data] 
    # indicator_all_zero_row = [1 not in row for row in matrix_without_missing_data] 
    #
    # # defining a local function
    # def local_fn_row_diseq(j):
    #   nr = indicator_all_zero_row[j] + indicator_all_one_row[j] - product_of_fractions_zeros - product_of_fractions_ones
    #   dr = 1 - product_of_fractions_ones - product_of_fractions_zeros
    #   return nr/dr
    # row_diseq = [local_fn_row_diseq(j) for j in range(0, n_rows)]
    #
    # END block comments

    # PartitionDisEq's input is already partitioned, whereas
    # HaplotypePartitionScore takes as input an unpartitioned haplotype,
    # and also a takes in a vector specifying the partition. In order to
    # use that HaplotypePartitionScore, I create the dummy_partition
    # vector.
    dummy_partition = [1 for k in range(0, n_cols)]
    row_diseq = [HaplotypePartitionScore(matrix_without_missing_data[j], number_of_ones, n_rows, dummy_partition) for j in range(0, n_rows)]
    matrix_diseq = sum(row_diseq)/n_rows

    return matrix_diseq

############################
############################

def TwoLocusLD(haplotype_matrix):
    """Compute LDs for each pair of sites. Thus, there are n_col choose 2
    of them, if n_col is the # columns in haplotype_matrix.

    For each pair of loci, discard data from an individual if there is
    missing data. Thus, different pairs of loci may be from different number
    of individuals.

    LD for a pair of sites i, j = f_{11} - p_i * p_j, where f_11 is the
    freq of the rows with 1 at both sites i and j. p_i and p_j are
    frequencies of type 1 at sites i and j respectively.
    """


    #sum over rows for each of the columns to get the derived allele
    #frequency in each column
    #derived alleles are denoted 1 in the hap matrix. 
    #Summing works because haplotypes are represented as 0-1 vectors.

    #print "HAPLOMAT =", haplotype_matrix
    (n_rows, n_col) = haplotype_matrix.shape


    #n_col = 10 PURELY FOR DEBUGGING PURPOSES
    #n_col = 10
    n_col_choose_2 = n_col*(n_col-1)/2

    #E is a list with n_col_choose_2 items. E[i] will the LD for the 2 loci 
    #extracted out by the mask represented set inside the while loop
    E = []

    k = n_col_choose_2-1
    bit_position_1 = 0
    bit_position_2 = 1

    while k > -1:
        print "analyzing partition", n_col_choose_2-k, "of ", n_col_choose_2
        partitioned_matrix = haplotype_matrix[0:n_rows, (bit_position_1, bit_position_2)]
        #print partitioned_matrix
        z = PartitionDisEq(partitioned_matrix, (bit_position_1, bit_position_2))
        E.append(z)
        k = k-1
        if (bit_position_2 == n_col-1):
            bit_position_1 = bit_position_1+1
            bit_position_2 = bit_position_1+1
        else:
            bit_position_2 = bit_position_2+1
          
    #reverse - because the lists are populated in decreasing order
    #according to the mask.          
    E.reverse()
    return E

############################
############################
