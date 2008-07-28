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

def list2string(l, replacement_dictionary = {}):
    """convert a list to a string, but convert any item that is a
       missingdata_token to a missingdata_repl
       """
    s = ''
    for item in l:
        if item in replacement_dictionary:
            item = replacement_dictionary[item]
        s = s + ' ' +str(item)
    s.strip()
    return s

############################
############################

def Standardize(h, missingdata_dict = {'#$%^': -1}):

    
#    print h, h.dtype
#    print "Dictionary = ", missingdata_dict
    for key in missingdata_dict:
        std_missingdata_token = missingdata_dict[key]


    (n_rows, n_col) = h.shape
    h.shape = n_rows*n_col
    for i in range(0, n_rows*n_col):
        if h[i] in missingdata_dict:
#            print "H, DIC[H] =", h[i], missingdata_dict[h[i]] 
            h[i] = missingdata_dict[h[i]]
    h.shape = (n_rows, n_col)            

    std_h = zeros(h.shape)

    biallelic_columns = []
    non_biallelic_columns = []
    for r in range(0, n_col):
        col = h[0:n_rows, r]
        #print col
        unique_items = list(set(col))
        
        #print "UNIQUE = ", r, unique_items
        #print col
        
        #if (std_missingdata_token in set(col)):
        #    assert len(unique_items) == 3, "Data not bi-allelic?" 
        #else:
        #    assert len(unique_items) == 2, "Data not bi-allelic?" 

        if (std_missingdata_token in set(col)):
            if len(unique_items) != 3:
                non_biallelic_columns.append(r)
                continue
        else:
            if len(unique_items) != 2:
                non_biallelic_columns.append(r)
                continue
        biallelic_columns.append(r)

        dic = {}
        dic[std_missingdata_token] = std_missingdata_token
        k = 0
        for i in range(0, len(unique_items)):
            if unique_items[i] == std_missingdata_token:
                continue
            dic[unique_items[i]] = k
            k = k+1

        std_col = [dic[item] for item in col]            
        std_h[0:n_rows, r] = std_col

    # At this point the non-biallelic columns, if any, would be all zeros. We now
    # filter out any such rows.
    if (len(non_biallelic_columns) > 0):
        print "Non-SNP columns: ", non_biallelic_columns
        std_h = std_h[0:n_rows, biallelic_columns]
    return std_h
                    

############################
############################

#def ConvertFromEIF(eif_mat):
#    """convert from the format of the EIF data that Yi-Ju provided. 
#       discard first five columns ( it is the id); then convert each row to two rows by
#       picking out alternate columns. 
#       """
#    
#    (n_rows, n_col) = eif_mat.shape
#    h = []
#
#    for i in range(0, n_rows):
#        row = eif_mat[i]
#        id = row[0:5]
#        genotype = row[5:n_col]
#
#        # even = [0, 2, 4, ..., n_col-5-2]
#        even = 2*arange((n_col-5)/2)
#        # odd = [1, 3, ..., n_col-5-1]
#        odd = 2*arange((n_col-5)/2)+1
#
#        #extract the even and odd elements from the row genotype
#        even_type = genotype[even]
#        odd_type = genotype[odd]
#
#        h.append(even_type)
#        h.append(odd_type)
#
#    h = array(h)
#    (x, y) = h.shape
#    assert (2*n_rows == x) and ((n_col-5)/2 == y)
#    return h

############################
############################

def Heterozygosity(genotype_mat, missingdata_token='#@$%'):
    """Take a matrix with phased or unphased genotype information; 
    that is the first two
    rows is the  genotype of the first individual, the next two rows
    the  genotype of the next individual and so on. This function
    creates a quote haplotype unquote matrix where the character at a locus
    denotes if the individual is heterozygous or homozygous at that locus. 

    Input: is a (2k x n) matrix, where k is the #individuals, and n is the
    #loci. Each pair of rows (2i, 2i+1), 0 < i < k denotes an
    unphased/phased genotype. Allelic states are denoted by single
    characters. Any characters can be used as long as there are exactly two
    characters per column. Different characters may be used for different
    columns; In other words the matrix need not be standardized. If one of
    the 2 bits is missing at a locus, then entire locus is considered
    missing.

    Output: A 0-1,-1 (assuming missingdata_token is -1) 
    matrix of dimension (k, n); the (i, j)th entry denotes if
    the i-th individual is heterozygous at site j or not. 1 denotes
    heterozygosity, and 0 denotes homozygosity.

    missingdata_token is set to a default of '@#$%', so that if it is not
    specified, the token will not match any entry in the input matrix.
    """
    
    n_rows = genotype_mat.shape[0]
    n_col =  genotype_mat.shape[1]

    assert n_rows%2 == 0, "Input Matrix Not Genotypes? There seems to be an odd \
    number of rows."

    trans = genotype_mat.transpose()

    #this is the tranpose of the matrix that will be returned.
    zygosity_matrix_transpose = zeros((n_col, n_rows/2))

    #note that n_col is the number of rows in the trasnposed matrix
    for r in range(0, n_col):
        col = list(trans[r])
        zygosity_col = zeros(n_rows/2)
        #zygosity_col = array([int(col[index] != col[index+1]) for index in range(0, n_rows, 2)])
        for index in range(0, n_rows, 2):
            if -1 in [col[index], col[index+1]]:
                zygosity_col[index/2] = -1
            else:
                zygosity_col[index/2] = int(col[index] != col[index+1])

        zygosity_matrix_transpose[r] = zygosity_col.copy()

    return zygosity_matrix_transpose.transpose()

############################
############################

def AlleleOneHomozygosity(genotype_matrix):
    """
    Takes a 2*M x N matrix G of genotypes and constructs a M x N character
    matrix Z, such that Z[i, j] = 1 if and only if the i^th zygote (rows
    2*i and 2*i+1) in G is homozygous for allele 1 at site j: i.e., Z[i, j]
    = 1 iff G[2*i, j] = G[2*i+1, j] = 1.

    Input parameters
    ----------------

    genotype_matrix     A (mostly) 0-1 numpy matrix 
                        (i.e., a numpy array of numpy arrays)
                        with dimensions of the form 
                        2*m x n, where m is interpreted as 
                        the number of zygotes, and n the number of sites.
                        
                        Each pair of rows (2i, 2i+1) is interpreted to denote 
                        a phased or unphased genotype. 

                        Apart from 0s and 1s, the matrix may also contain
                        -1 denoting missing data. 

    Return value
    ------------

    allele_one_homozygosity     A numpy m x n matrix, such that 
                                
                                allele_one_homozygosity[i, j] = 1 if and
                                only if genotype_matrix[2*i, j] =
                                genotype_matrix[2*i+1, j] = 1 

                                allele_one_homozygosity[i, j] = -1 if 
                                genotype_matrix[2*i, j] = -1 or
                                genotype_matrix[2*i+1, j] = -1.

                                allele_one_homozygosity[i, j] = 0 in all
                                other cases.
    """

    (n_rows, n_cols) = genotype_matrix.shape

    assert n_rows%2 == 0, "Input Matrix Not Genotypes? There seems to be an odd \
    number of rows."
    
    n_genotypes = n_rows/2

    allele_one_homozygosity = zeros(n_genotypes, n_cols)
    
    for j in range(0, n_genotypes):
        # kind of obscurantist coding. But I'll explain. 
        #
        # 'minimum' is a numpy function that makes element-by-element
        # comparison of two vectors and puts the minimum of each comparsion
        # into a vector. That is:
        # allele_one_homozygosity[i] = min(genotype_matrix[2*j, i], current_zygote[2*j+1, i])
        #
        # Now, this is what we want allele_one_homozygosity to be, since
        # min(1, 1) = 1, min(0, 1) = 0, and min(-1, 1) = min(-1, 0) = min(-1, -1) = -1.
        allele_one_homozygosity[j] = minimum(genotype_matrix[2*j], genotype_matrix[2*j+1])
    
    return allele_one_homozygosity


############################
############################

def FilterIncompleteRecords(haplotype_matrix, missingdata_token):
    """Filter out any row with missing data.
    """

    (n_rows, n_col) = haplotype_matrix.shape
    filtered_matrix = []
    for j in range(0, n_rows):
        if missingdata_token in haplotype_matrix[j]:
            continue
        filtered_matrix.append(list(haplotype_matrix[j]))

    return array(filtered_matrix)

############################
############################
