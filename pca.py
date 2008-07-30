#!/usr/bin/env python

# Using CapWords for function, method, and class names
# Using underscored_names for variable names
# Using firstWordLower for module, package names.
# Using ALL_CAPS_WITH_UNDERSCORES for file handles

#standard imports 
import sys

#python extensions
import numpy 

#Ganesh's modules
import readData 
import formatConversions 
import ldComputations
import levelDistancePartition
import fileWrite

"""
The script pca.py takes a haplotype matrix H (dimensions M x N), 
and outputs a matrix of disequilibrium measures, V, of dimension M x (2^N - 1). 

a. All disequilibrium measures computed are single-haplotype versions of
  the multilocus disequilibrium measure E (see bottom of page 2 of
  Dpca.pdf). 

  For a particular haplotype, h, and a set of sites Z, E(h, Z) is calculated as:

  set: f_1 = 1 if h[j] = 1 in all sites j in Z
           = 0 otherwise
       
       f_0 = 1 if h[j] = O in all sites i in Z
           = 0 otherwise.

       P_1(Z) = product of allele-type 1 sample frequencies over all
                sites in Z. Frequencies here means fractions, not
                actual numbers. 

       P_0(Z) = product of allele-type 0 sample frequencies over all
                sites in Z. 

       Then,

                   f_1 + f_0  - P_1(Z) - P_0(Z)
       E(h, Z) =   ----------------------------
                      1 - P_1(Z) - P_0(Z)


b. V[i] (i.e., row i in V) is a vector of disequilibrium measures calculated for the
  haplotype H[i] (i.e., row i in H). 
  
  V[i, j] (i.e., the j-th element of
  V[i]) is the multilocus disequilibrium computed for the partition of
  H[i] represented by the binary representation of j+1. 

  For example, let N = 4 and j = 5.
  Then, V[i, 5] is the multilocus disequilibrium of the
  partition of haplotype specified by [0, 1, 1, 0] i.e., of the partition
  [H[i, 1], H[i, 2]. So in this example, V[i, 5] is the traditional 2-locus linkage
  disequilibrium between sites 1 & 2, computed for the *single haplotype H[i]*

c. pca.py should be invoked as:
   "python pca.py <parameters>", where the parameters are as follows:
                                        
   Parameters (parameters in [] optional)
   --------------------------------------
   -help
   -hapfile <haplotype matrix file>
   -missingdata <missingdata token> 
   -output <output file>
   -random <a random integer j such that 0 <= j < number of sites in input matrix>

   For example: 
   1. python pca.py -hapfile nomix1214.pop1.seq -output vectors_for_pca.test1.out

       reads the haplotype matrix from the file nomix1214.pop1.seq and writes the
       output in the file vectors_for_pca.test1.out

       If a token other than -1 is used to denote missing data, then it
       must be specified after the '-missingdata' argument. By default,
       -1 is assumed to denote missing data.

   2. python pca.py -hapfile nomix1214.pop1.seq -output vectors_for_pca.test1.out -random 10

        similar to one, but randomly samples 10 sites, and uses the sampled
        sites as the haplotype matrix. In the new matrix, the sites appear
        in the same order in which they appeared in the original matrix.
        It also outputs which sites were sampled. 
        
"""
       

def IsOption(s):
    """does string s start with a - ? like -Hap"""
    return s[0] == '-'


OP_FILE = None
random_sites = False
n_provided_options = 0
non_std_missingdata_token = None
std_missingdata_token = '-1'

arguments = levelDistancePartition.PadNone(sys.argv)

# discard the first argument; this argument is this source python file's name
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

    if (arg == '-output'):
        output_file_name = arguments.next()
        # open for append, but clear out the buffer first.
        OP_FILE = open(output_file_name, 'w')
        OP_FILE.close()
        OP_FILE = open(output_file_name, 'a')
        arg = arguments.next()
        continue

    if (arg == '-random'):
        random_sites = True
        number_of_random_sites = int(arguments.next())
        arg = arguments.next()
        continue

    if (arg == '-help'):
        arg = arguments.next()
        usage = """Usage: python pca.py <parameters>
                                            
                                               Parameters (parameters in [] optional)
                                               --------------------------------------
                                               -help
                                               -hapfile < haplotype matrix file >
                                               -missingdata < missingdata token > 
                                               -output < output file>
                                               -random <a random integer j such that 0 <= j < number of sites in input matrix>
                                              """
        print usage
        sys.exit()

# read the data matrix from a file.
# non_std_mat is  2-D numpy array, but the datatype is *NOT INT*
nonstd_haplotype_matrix = readData.ReadFromTxt(haplotype_file_name)

# make 0 and 1 the two allelic types; after replacing the
# missingdata_token first. 
if not non_std_missingdata_token == None:
    print 'Missing data token = ', non_std_missingdata_token
    standard_haplotype_matrix = formatConversions.Standardize(nonstd_haplotype_matrix, dict([(non_std_missingdata_token, std_missingdata_token)]))
    
else:
    print 'Missing data token ', 'not provided: assuming there is no missing data'
    standard_haplotype_matrix = formatConversions.Standardize(nonstd_haplotype_matrix)

(number_of_haplotypes, number_of_sites) = standard_haplotype_matrix.shape

if random_sites == True:
    list_of_sites = []
    for j in range(0, number_of_random_sites):
        # Return random integers x such that 0 <= x < number_of_sites.
        random_number = numpy.random.randint(0, number_of_sites)
        while random_number in list_of_sites:
            random_number = numpy.random.randint(0, number_of_sites)
        list_of_sites.append(random_number)

    list_of_sites.sort()        
    if OP_FILE == None:
        print "List of sites picked (site numbers start from 0)"
        print list_of_sites
        print "\n"
    else:
        OP_FILE.write("List of sites picked (site numbers start from 0)\n")
        OP_FILE.write(str(list_of_sites))
        OP_FILE.write("\n\n")
    standard_haplotype_matrix = standard_haplotype_matrix[0:number_of_haplotypes, list_of_sites]
    

no_missing_data_standard_haplotype_matrix = ldComputations.FilterMissingData(standard_haplotype_matrix)
number_of_rows_no_missing_data = (no_missing_data_standard_haplotype_matrix.shape)[0]
allele_one_site_frequencies = ldComputations.AlleleOneSiteFrequencies(no_missing_data_standard_haplotype_matrix)

# for testing and debugging purposes:
#no_missing_data_standard_haplotype_matrix = no_missing_data_standard_haplotype_matrix[0:number_of_rows_no_missing_data, 0:4]

# vectors_for_pca[j] =
# ldComputations.PerHaplotypeScoresVector(no_missing_data_standard_haplotype_matrix[j], 
#                                         allele_one_site_frequencies, 
#                                         number_of_rows_no_missing_data)

vectors_for_pca = [ldComputations.PerHaplotypeScoresVector(no_missing_data_standard_haplotype_matrix[j], allele_one_site_frequencies, number_of_rows_no_missing_data) for j in range(0, number_of_rows_no_missing_data)]

if OP_FILE == None:
    print "Vectors for PCA. One line per haplotype\n"
    print vectors_for_pca
else:
    OP_FILE.write("Vectors for PCA. One line per haplotype\n\n")
    fileWrite.Print2DArrayToFile(OP_FILE, numpy.array(vectors_for_pca))
