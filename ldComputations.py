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

############################
############################

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
    row_diseq.sum()/Q) has the advantage that it can be used without
    modification to calculate Marcy's individual haplotype scores (1. See
    Section 2, Dpca.pdf. 2. Marcy also defines similar individual zygote
    scores, but PartitionDisEq cannot be used to calculate the zygote
    scores, except the ID disequilibrium.)
    """

    (n_rows, n_col) = partition_matrix.shape

    if n_col < 2:
        print 'Cannot calculate disequilibrium for the single column ', partition
        return -1

    # disregard rows with -1 (missing data)
    matrix_without_missing_data = filter(lambda row: -1 not in row, partition_matrix)
    # the statement above returns a normal list of numpy arrays. The following
    # statement typecasts the normal list into a numpy array,
    # like partition_matrix is. Note that we are not modifying
    # partition_matrix here at all.
    matrix_without_missing_data = array(matrix_without_missing_data)

    (n_rows, n_cols) = matrix_without_missing_data.shape

    # each row in matrix_without_missing_data.transpose() is a column n
    # matrix_without_missing_data. Also note: calling the transpose
    # method of the array object does not modify the object itself.
    # also, we would like number_of_ones and number_of_zeros to be numpy
    # arrays.
    number_of_ones = array([column.sum() for column in matrix_without_missing_data.transpose()])
    # n_rows - number_of_ones = [n_rows-x for x in number_of_ones]
    number_of_zeros = n_rows - number_of_ones

    product_of_fractions_ones  = float(number_of_ones.prod())/(n_rows**n_cols)
    product_of_fractions_zeros = float(number_of_zeros.prod())/(n_rows**n_cols)

    indicator_all_one_row = [0 not in row for row in matrix_without_missing_data] 
    indicator_all_zero_row = [1 not in row for row in matrix_without_missing_data] 

    # defining a local function
    def local_fn_row_diseq(j):
        nr = indicator_all_zero_row[j] + indicator_all_one_row[j] - product_of_fractions_zeros - product_of_fractions_ones
        dr = 1 - product_of_fractions_ones - product_of_fractions_zeros
        return nr/dr
    row_diseq = [local_fn_row_diseq(j) for j in range(0, n_rows)]
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

#def MultiLocusLD(haplotype_matrix):
#    """H = set of loci; for each subset B of $2^H$ except the emptyset,
#       calculate $f_{B, 1} and $f_{B, 0}$. $f_{B, 1}$ is the frequency of
#       haplotypes that have all $1$s in the loci in the subset $B$. Likewise
#       for $f_{B, 0}$. Then E_{B} for the set of loci B is calculated as 
#       E_{B} = \frac{f_{B, 1} - \Pi_{j \in B} p_j + f_{B, 0} - \Pi_{j \in B} (1-p_j)}{1 - \Pi_{j \in B} p_j - \Pi_{j \in B} (1-p_j)}
#
#       General algorithm: say for $n (#genes) = 4$, apply masks from $0001$ to 
#       $1111$ to the columns. Each mask represents a subset $B$. 
#       The columns that are 1 in the mask are extracted,
#       logically speaking; Then count the #rows whose masked sum = #ones in the
#       mask. This would be $f_{B, 1}$. Similarly, calculate $f_{B, 0}. Also
#       calculate p_j for each column in the partition. All this
#       is implemented in PartitionDisEq.
#
#       Input: haplotype_matrix, an (k, n) 0-1 numpy array; each row represents the
#       haplotype of an individual.
#
#       Output: a list of 2**n - 1 LD measures,  one for each
#       non-empty subset of the n loci.
#       """
#
#    (n_rows, n_col) = haplotype_matrix.shape
##
#    #E is a list with 2**n-1 items. E[i] will the LD for the loci extracted
#    #out by the mask represented by i.
#    E = []
#
#    #n_col = 10 PURELY FOR DEBUGGING PURPOSES
#    #n_col = 10
#    k = 1
#    two_to_the_n_col = 2**n_col
#    # had to use  while loop instead of for because the range function
#    # failed with large numbers such as 2**n_col
#    while k < two_to_the_n_col:
#        print "analyzing partition", k, "of ", two_to_the_n_col
#        #BitArray returns a numpy array
#        mask = BitArray(k, n_col)
#
#        #mask_sum = # ones in the mask.
#        mask_sum = mask.sum()
#        if (mask_sum != 2):
#            k = k+1
#            continue
#
#        #partition contains the positions in the mask that are set to 1
#        w = where(mask == 1)
#        partition = w[0]
#        partitioned_matrix = haplotype_matrix[0:n_rows, partition]
#        z = PartitionDisEq(partitioned_matrix, partition)
#        E.append((w, z))
#        #E.append(z)
#        k = k+1
#            
#    return E
#
#############################
###############################
