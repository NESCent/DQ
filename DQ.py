#!/usr/bin/env python

# Using CapWords for function, method, and class names
# Using underscored_names for variable names
# Using firstWordLower for module names;
# Using CapWords for package names;
# Using ALL_CAPS_WITH_UNDERSCORES for file handles

#standard imports 
import sys
import os
import re

#python extensions
from numpy import * 

#Ganesh's modules
from readData import *
from formatConversions import *
from levelDistancePartition import *
from fileWrite import *

def IsOption(s):
    """does string s start with a - ? like -Hap"""
    return s[0] == '-'


#this option *must* be provided by the user
haplotype_file_name = None

#default values for optional parameters
idd_ld = 'ld'
is_allpairs = 'no'
OP_FILE = sys.stdout
output_file_name = OP_FILE.name

#these options are optional, but are set to None so that it can be tested
#if they are user-specified.
level = fixed_distances = partition = None
non_std_missingdata_token = None

n_provided_options = 0
std_missingdata_token = '-1'

arguments = PadNone(sys.argv)

# eat the first argument; this argument is this source python file's name
arg = arguments.next()
arg = arguments.next()
while not(arg == None):
    assert IsOption(arg), 'Unknown parameter '+str(arg)

    if (arg == '-hapfile'):
        haplotype_file_name = arguments.next()
        arg = arguments.next()
        continue

    if (arg == '-missingdata'):
        non_std_missingdata_token = arguments.next()
        arg = arguments.next()
        continue

    if (arg == '-allpairs'):
        is_allpairs = arguments.next()
        arg = arguments.next()
        continue

    if (arg == '-idd_ld'):
        idd_ld = arguments.next()
        arg = arguments.next()
        continue

    if (arg == '-output'):
        output_file_name = arguments.next()
        # open for append, but clear out the buffer first.
        OP_FILE = open(output_file_name, 'w')
        OP_FILE.close()
        OP_FILE = open(output_file_name, 'a')
        arg = arguments.next()
        continue

    if (arg == '-level'):
        level = int(arguments.next())
        n_provided_options += 1
        arg = arguments.next()
        continue

    if (arg == '-dist'):
        fixed_distances = []
        n_provided_options += 1
        arg = arguments.next()
        while not(arg == None or IsOption(arg)):
            fixed_distances.append(int(arg))
            arg = arguments.next()
        continue

    if (arg == '-partition'):
        partition = []
        n_provided_options += 1
        arg = arguments.next()
        while not(arg == None or IsOption(arg)):
            partition.append(int(arg))
            arg = arguments.next()
        continue

    if (arg == '-help'):
        arg = arguments.next()
        usage = """Usage: python DQ.py < parameters >
                                            
                                               Parameters (parameters in [] optional)
                                               --------------------------------------
                                               -help
                                               -hapfile < haplotype matrix file >
                                               -missingdata < missingdata token > -output < output file > 
                                               -allpairs <yes/no> 
                                               -idd_ld < idd or ld >
                                               -output < output file>
                                              [-level < number of loci >] 
                                              [-dist < the list of fixed distance to be maintained between loci >] 
                                              [-partition < a list of loci
                                              for which the disequilibrium has to be calculuated. Assumed that loci are 
                                              numbered from 0 through n_col-1 >] 
                                              """
        print usage
        sys.exit()


assert n_provided_options < 2, 'Too many options: ' + usage

print 'Input File = ', haplotype_file_name
print 'Output file = ', output_file_name
print 'IDD or LD? ', idd_ld
print 'Two-locus disequilibria in matrix form?', is_allpairs


# read the data matrix from a file.
# non_std_mat is  2-D numpy array, but the datatype is *NOT INT*
nonstd_mat = ReadFromTxt(haplotype_file_name)

# make 0 and 1 the two allelic types; after replacing the
# missingdata_token first. 
if not non_std_missingdata_token == None:
    print 'Missing data token = ', non_std_missingdata_token
    haplotype_matrix = Standardize(nonstd_mat, dict([(non_std_missingdata_token, std_missingdata_token)]))
else:
    print 'Missing data token = ', 'not provided: assuming there is no missing data'
    #haplotype_matrix = Standardize(nonstd_mat, dict([(non_std_missingdata_token, std_missingdata_token)]))
    haplotype_matrix = Standardize(nonstd_mat)

if (idd_ld == 'idd'): 
    data_matrix = Heterozygosity(haplotype_matrix)
else:
    data_matrix = haplotype_matrix

tolerance = 0
(n_rows, n_col) = data_matrix.shape


if n_provided_options == 0 or is_allpairs == 'yes':
    # by default compute all pairwise disequilibria and write them as a
    # matrix
    E_Z = TwoLocusLD(data_matrix)
    #n_col = 10
    e_z_2locus = Make2dArray(E_Z, n_col)
    #OP_FILE.write('**************OUTPUT***************\n')
    Print2LocusLDArrayToFile(OP_FILE, e_z_2locus)

elif not(level == None):
    print 'chosen level (= chosen #loci), Total #loci = ', level, n_col
    # level = #loci; for level = 2, pairwise diseq are computed; but not
    # all pairwise diseq. only a log(n_col) # distances: 1, 2, 4, ..., n/8,
    # n/4, n/2 are used. 
    #for now level can only be 2
    assert (level in range(1, n_col+1)), 'Illegal level: level should be between 1 and ' + str(n_col)
    [distances, disequilibria] = AllFixedDistances(data_matrix, level, tolerance)
    #OP_FILE.write('**************OUTPUT***************\n')
    PrintFixedDistDiseqToFile(OP_FILE, distances, disequilibria)

elif not(fixed_distances == None):
    print 'inter-locus distances, total #loci', fixed_distances, n_col
    # the fixed_distances could be basically a set of distances. If there
    # is only one distance in the set, then pairwise diseq are computed for
    # those sites that are the given distance apart. In general, the
    # distance could be (r1, r2, ... r_k); then (k+1)-way disequilibria
    # should be computed. For each set of (k+1) sites, the distance between 
    # the i-th and the (i+1)th site is r_i
    lds = OneFixedDist(data_matrix, fixed_distances)
    #OP_FILE.write('**************OUTPUT***************\n')
    OP_FILE.write(str(fixed_distances) + ' ' + str(lds) + '\n')
    #print fixed_distances, lds

elif not(partition == None):
    print 'Partition, total #loci', partition, n_col
    # a specific partition: i.e., a list of loci (numbered from 0 through
    # n_col-1), for which the disequilibrium has to be calculated.
    ld = OnePartition(data_matrix, partition)
    #OP_FILE.write('**************OUTPUT***************\n')
    OP_FILE.write(str(ld) + '\n')
    #print ld

OP_FILE.close()    
