# Copyright (C) 2013 Kyle Hernandez

# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program; if not, see <http://www.gnu.org/licenses/>.

import time
import sys
from itertools import izip, tee
import gzip

def main():
    """
    ---------------------------------------------------------------------------
    AUTHOR: Kyle Hernandez
    EMAIL: kmhernan@utexas.edu

    Calculate the distances between snps on each chromosome.
    ---------------------------------------------------------------------------

    USAGE: python distances_between_snps.py gmatrix.tab out.tab

    ARGUMENTS:
	gmatrix.tab - Genotype Matrix file
        out.tab     - Name of output file
    """
    snp_dict = load_gmatrix()
    get_distances(snp_dict)

def load_gmatrix():
    """
    Parses the gmatrix file into a dictionary
    """
    dic = {}
    try:
        with gzip.open(gmatrix, 'rb') as f:
            for line in f:
                if line.startswith('L'): continue
            
                cols  = line.rstrip().split('\t')
                scaff = cols[0].split(':')[0]
                pos   = int(cols[0].split(':')[1])

                if scaff not in dic: dic[scaff] = {}
                dic[scaff][pos] = ''
        return dic

    except:
        with open(gmatrix, 'rU') as f:
            for line in f:
                if line.startswith('L'): continue
            
                cols = line.rstrip().split('\t')
                scaff = cols[0].split(':')[0]
                pos   = int(cols[0].split(':')[1])

                if scaff not in dic: dic[scaff] = {}
                dic[scaff][pos] = ''
        return dic

def get_distances(snp_dict):
    """
    Calculates the distances between consecutive snps on a chromosome.
    """
    with open(outfil, 'wb') as o:
       o.write("Chromosome\tDistance\n")
       for c in sorted(snp_dict):
           curr_pos = sorted(snp_dict[c].keys())
           curr_diff = [y-x for x,y in pairwise(curr_pos)]
           for i in curr_diff: 
               o.write(c + '\t' + str(i) + '\n')

def pairwise(iterable):
    """
    Python cookbook pairwise comparison of iterable elements
    s -> (s0,s1), (s1,s2), (s2, s3),...
    """
    a, b = tee(iterable)
    next(b, None)
    return izip(a, b) 

if __name__=='__main__':
    start = time.time()
    if len(sys.argv) != 3:
        print main.__doc__
        sys.exit()
    gmatrix = sys.argv[1]
    outfil  = sys.argv[2]
    main()
    print "Finished; Took:", time.time() - start, "seconds."
