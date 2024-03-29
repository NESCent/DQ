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

#Ganesh's modules
from formatConversions import *
from readData import *

############################
############################

def PrintAllPartitionsDisEq(OP_FILE, E_Z):
    for item in E_Z:
        OP_FILE.write(str(item))
        OP_FILE.write('\n')

def Print2DArrayToFile(OP_FILE, m):
    
    for i in range(0, m.shape[0]):
        for j in range(0, m.shape[1]):
            #OP_FILE.write('%9.7f' % m[i, j])
            OP_FILE.write(str(m[i, j]))
            OP_FILE.write('\t')
        OP_FILE.write('\n')

def Print2LocusLDArrayToFile(OP_FILE, m):
    
    OP_FILE.write('\t')
    for i in range(0, m.shape[0]):
        s = 'SNP'+str(i+1)+'\t'
        OP_FILE.write(s)
    OP_FILE.write('\n')
    for i in range(0, m.shape[0]):
        s = 'SNP'+str(i+1)+'\t'
        OP_FILE.write(s)
        for j in range(0, m.shape[1]):
            OP_FILE.write('%9.7f' % m[i, j])
            #OP_FILE.write(str(m[i, j]))
            OP_FILE.write('\t')
        OP_FILE.write('\n')


############################
############################

def Make2dArray(e_list, n_bits):
    """Takes a list of pairwise LD values)
    and puts them in the form of an array e_2d_array. The size of e_list
    must be (n_bits choose 2). This function assumes that e_list is
    generated by TwoLocusLD. e_list[0] corresponds to the smallest mask
    with two ones. E.g., for n_bits = 5, e_list[0] corresponds to the mask
    00011; therefore e_list(0) has the LD between loci 3 & 4 (count starts
    from 0. The left most bit is the 0th locus). 

    Mask in decimal and binary:  e_list index
    3 [ 0.  0.  0.  1.  1.]         0
    5 [ 0.  0.  1.  0.  1.]         1
    6 [ 0.  0.  1.  1.  0.]         2
    9 [ 0.  1.  0.  0.  1.]         3
    10 [ 0.  1.  0.  1.  0.]        4
    12 [ 0.  1.  1.  0.  0.]        5
    17 [ 1.  0.  0.  0.  1.]        6
    18 [ 1.  0.  0.  1.  0.]        7
    20 [ 1.  0.  1.  0.  0.]        8
    24 [ 1.  1.  0.  0.  0.]        9

    Note that going through e_list backwards corresponds to a logical way
    of cycling through the bits of the mask. This is what allows us to
    infer the mask even though the mask is not part of the input.
    """

    n_bits_choose_2 = n_bits*(n_bits-1)/2
    ld_2locus = zeros((n_bits, n_bits))
    j = n_bits_choose_2-1
    bit_position_1 = 0
    bit_position_2 = 1
    while j > -1:
        ld_2locus[bit_position_1, bit_position_2] = e_list[j]
        j = j-1
        if (bit_position_2 == n_bits-1):
            bit_position_1 = bit_position_1+1
            bit_position_2 = bit_position_1+1
        else:
            bit_position_2 = bit_position_2+1

    return ld_2locus            

def PrintFixedDistDiseqToFile(OP_FILE, distances, disequilibria):

    assert (len(distances) == len(disequilibria)), 'Unequal # distances and disequilibria'

    for d, e in zip(distances, disequilibria):
        s = str(d) + ' ' + str(e) + '\n'
        OP_FILE.write(s)


