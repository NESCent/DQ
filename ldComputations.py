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

    (n_rows, n_col) = partition_matrix.shape

    if n_col < 2:
        print 'Cannot calculate disequilibrium for the single column ', partition
        return -1

    freq0 = 0.0
    freq1 = 0.0
    freq_check = 0.0
    n_complete_rows = 0
    
    #one for each locus 
    derived_allele_freq = zeros(n_col) 
    #one for each locus 
    anc_allele_freq = zeros(n_col)

    for r in range(0, n_rows):
        row = partition_matrix[r]
        if -1 in row:
            continue

        n_complete_rows += 1

        #print row
        s = row.sum()

        #numpy allows such cool matrix additions. it even understands
        #1-row!

        derived_allele_freq += row
        #1-row is row bit flipped!
        anc_allele_freq += 1-row

        if s == n_col:
            #all ones
            freq1 += 1
        elif s == 0:
            #all zeros
            freq0 += 1
        else:
            #check sum
            freq_check += 1

    derived_allele_sum = derived_allele_freq.sum()
    derived_allele_prod = derived_allele_freq.prod()
    anc_allele_prod = anc_allele_freq.prod()

    assert (freq0+freq_check+freq1 == n_complete_rows)

    #The following two asserts must hold *ONLY FOR* the 2-Locus case. That
    #is, *ONLY WHEN PARTITION_MATRIX HAS 2 Columns*

    #assert (freq_check+2*freq1 == derived_allele_sum), "Frequencies don't add up!"
    #assert (n_complete_rows*freq1)-derived_allele_prod == (n_complete_rows*freq0)-anc_allele_prod, 'Mismatch between two ways of calculating LD'

    assert (freq0 + freq1 <= n_complete_rows), 'Frequencies add up to more than 1!'

    # The following three lines are in error. 
    #   x = (n_complete_rows*freq1-derived_allele_prod)
    #   y = (n_complete_rows*freq0-anc_allele_prod)
    #   u = n_complete_rows**2-derived_allele_prod-anc_allele_prod
    # end erroneous code

    # The erroneous lines are being replaced with the following lines
    x = ((n_complete_rows**(n_col-1))*freq1-derived_allele_prod)
    y = ((n_complete_rows**(n_col-1))*freq0-anc_allele_prod)
    u = (n_complete_rows**n_col)-derived_allele_prod-anc_allele_prod
    # end replacement code

    if (x+y == u == 0):
        z = 1
        print "Something fishy among columns ", partition, ". Probably the \
        columns are all ancestral or all descendant allele types?"
    else:
        z = (x+y)/u

    #print n_complete_rows, freq0, freq1, derived_allele_prod, anc_allele_prod
    #print x, y, x+y, u
    assert (z <= 1 or (x+y < 0 and u < 0)) 

    #print '-------------'
    #print n_complete_rows, derived_allele_freq, anc_allele_freq, derived_allele_prod, anc_allele_prod 
    #print freq0, freq1, freq_check 
    #print 'x = ', x, '\t y = ', y, '\t u = ', u, '\t z = ', z 
    #print '-------------'

    return z
        

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
