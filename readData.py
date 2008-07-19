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

def max(x, y):
    if x > y:
        return x
    else:
        return y

############################
############################
def ReadFromTxt(txt_file):

    print 'TXT = ', txt_file
    TXT_FILE = open(txt_file, "r")
    #calculate the #rows (#individuals) and the #columns (#genes) in the
    #sample.
    n_rows = 0
    h = []
    max_length = 0
    for line in TXT_FILE.readlines():
        x = line.split()
        #print x, n_rows
        lengths = [len(s) for s in x]
        m = reduce(max, lengths)
        if (m > max_length):
            max_length = m
        #convert into integers.
        #int_x = [int(i) for i in x]
        h.append(x)
        n_col = len(x)
        n_rows = n_rows+1
    TXT_FILE.close()

    if max_length == 1:
        h = array(h, dtype = 'S2')
    else:
        h = array(h)
    return h

############################
############################

#def ReadFromFastphaseOutput(fastphase_file):
#    """Read input in fastphase format and return just the data matrix as as
#       numpy array
#
#       skip lines until 'BEGIN GENOTYPES'
#       thenm skip lines that start with #  These are id tags.
#       in the input each gene would have two states, but they can be denoted by
#       anything. 
#
#       Input:   the name of the data file (a fastphase output file)
#       Output:  a 2-d numpy array. datatype '|S1' (one-charcater string?)
#       data in the input file. *NOTE THAT THE DATATYPE IS NOT INT OR BOOL*
#       """
#
#    FF_FILE = open(fastphase_file, "r")
#
#    n_rows = 0
#    nonstd_mat = []
#    flag = False
#    hash = r'^#'
#    re_hash = re.compile(hash)
#
#    for line in FF_FILE.readlines():
#        if (line.strip() == 'BEGIN GENOTYPES'):
#            flag = True
#            continue
#        if (flag == True):
#            if re_hash.match(line) != None or len(line.strip()) == 0  or line.strip() == 'END GENOTYPES':
#                continue
#            else:
#                n_col = len(line.split()) 
#                n_rows = n_rows+1
#                nonstd_mat.append(line.split())
#    FF_FILE.close()
#    return array(nonstd_mat)                

############################
############################
