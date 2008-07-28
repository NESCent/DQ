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

   For example: 
   python pca.py -hapfile nomix1214.pop1.seq -output vectors_for_pca.test1.out

   reads the haplotype matrix from the file nomix1214.pop1.seq and writes the
   output in the file vectors_for_pca.test1.out

   If a token other than -1 is used to denote missing data, then it
   must be specified after the '-missingdata' argument. By default,
   -1 is assumed to denote missing data.
