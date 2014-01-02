#!/usr/bin/env python
# Filter sites with N alleles != 2; convert to 0/1.
# Copyright (C) 2014 copyright holder

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


import gzip, sys, time

def main():
    """
    ---------------------------------------------------------------------------
    AUTHOR: Kyle Hernandez
    EMAIL: kmhernan@utexas.edu

    Filter sites with n alleles != 2.
    Converts the character genotypes into 0 and 1 integers.
    ---------------------------------------------------------------------------
    USAGE: python get_biallelic_numeric.py gmatrix.tab ref.fa out.tab

    ARGUMENTS:
    	gmatrix.tab - Tab delimited genotype matrix file.
   	ref.fa      - Reference fasta file.
     	out.tab     - Name of output genotype matrix file.
    """
    site_dict, header = get_biallelic_sites()
    process_reference(site_dict, header)

def get_biallelic_sites():
    """
    Remove loci with n alleles != 2
    """
    dic = {}
    header = ''
    with gzip.open(gmatrix, 'rb') as f:
        for line in f:
            if line.startswith('L'):
               header = line

            else:
                cols = line.rstrip().split('\t')
                gts  = cols[1::]

                if len(set([i for i in gts if i != 'N'])) == 2:
                    scf, po = cols[0].split(':')
                    po = int(po) - 1

                    if scf not in dic: dic[scf] = {}
                    dic[scf][po] = gts
    return dic, header

def process_reference(dic, hd):
    """
    Find the genotype call in the reference and convert the genotypes
    to 0 if ref, 1 if not, 'N' if 'N'
    """
    scf = ''
    seq = []
    o = open(omatrix, 'wb')
    o.write(hd)
    with open(ref, 'rU') as f:
        for line in f:
            if line.startswith('>') and not scf:
                scf = line.rstrip()[1::]
            elif line.startswith('>') and scf:
                try:
                    curr_sites = dic[scf]
                    cseq = ''.join(seq)
                    # Use nested list comps to return a list of 
                    # tuples (position, numeric gts) 
                    conv = [(i, ['N' if j == 'N' else '0' if j == cseq[i] else
                                 '1' for j in curr_sites[i]])
                           for i in curr_sites.keys()]
                    print_conv(scf, conv, o)
 
                except KeyError: pass 

                scf = line.rstrip()[1::]
                seq = []
            else:
                seq.append(line.rstrip())
        try:
            curr_sites = dic[scf]
            cseq = ''.join(seq)
	    # Use nested list comps to return a list of 
            # tuples (position, numeric gts) 
            conv = [(i, ['N' if j == 'N' else 
      	                 '0' if j == cseq[i] else
                         '1' for j in curr_sites[i]])
                    for i in curr_sites.keys()]
            print_conv(scf, conv, o)
 
        except KeyError: pass 
    o.close()

def print_conv(scf, conv, o):
    """
    Prints the current scaffold data from the list of tuples to
    the omatrix stream.
    """
    [o.write('{0}:{1}\t{2}\n'.format(scf, i[0], '\t'.join(i[1])))
     for i in conv]

if __name__ == '__main__':
    start = time.time()
    if len(sys.argv) != 4:
        print main.__doc__
        sys.exit()
    gmatrix = sys.argv[1]
    ref     = sys.argv[2]
    omatrix = sys.argv[3]
    main()
    print 'Finished; Took: {: 0.5F} seconds.'.format(time.time() - start)
