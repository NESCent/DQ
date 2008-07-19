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


def ExtractNumbers(s):
    """given a string s = '[1, 2, 3]' (or similar), return the list ['1',
    '2', '4']"""

    t = s.strip('[]\n')
    comma_space = r', '
    re_comma_space = re.compile(comma_space)
    z = re_comma_space.split(t)
    #print z
    return z


dq_file = sys.argv[1] 
r_file = sys.argv[2]
print dq_file, r_file

DQ_FILE = open(dq_file, "r")
R_FILE = open(r_file, "w")
for line in DQ_FILE.readlines():
    x = line.split('] [', 1)
    print x
    coords = ExtractNumbers(x[0])
    print coords
    diseq_values = ExtractNumbers(x[1])
    print diseq_values
    n_diseq_values = len(diseq_values)
    for j in range(0, n_diseq_values):
        output_line = ""
        for k in range(0, len(coords)):
            output_line = output_line + str(coords[k]) + ' '
        output_line = output_line + str(diseq_values[j]) + '\n' 
        R_FILE.write(output_line)    
DQ_FILE.close()        
R_FILE.close()
