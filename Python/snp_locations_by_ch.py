#!/usr/bin/env python
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
from multiprocessing import Pool
import gzip
from collections import OrderedDict
from itertools import groupby
from operator import itemgetter

gff_dict  = OrderedDict() 
snp_dict  = OrderedDict() 

def main():
    """
    ---------------------------------------------------------------------------
    AUTHOR: Kyle Hernandez
    EMAIL: kmhernan@utexas.edu

    Calculate the distribution of polymorphic RAD loci across windows.
    ---------------------------------------------------------------------------

    USAGE: python snp_locations_by_window.py gmatrix.tab file.gff out.tab n_threads wsize

    ARGUMENTS:
    	gmatrix.tab - Tab-delimited genotype matrix file of variant sites
        file.gff    - GFF file
        out.tab     - Output file of counts
        n_threads   - The number of threads to run
        wize        - Size of window
    """

    # Load the GFF and SNP positions into dictionaries
    load_gff()
    process_matrix()
    
    # Map:
    # Create a pool of n_threads workers and use them to process
    # scaffolds separately
    ch_vals = gff_dict.keys()
    sys.stdout.write("Counting features...\n")
    pool    = Pool(processes=n_threads)
    ct_list = pool.map(process_dicts, ch_vals)

    # Reduce:
    # Process the list of dicts
    print_counts(ct_list)

def load_gff():
    """
    Parses GFF file into a dictionary.
    I should probably cpickle this (or JSON)
    """
    sys.stdout.write("Loading GFF dictionary...\n")
    with open(gff_fil, 'rU') as f:
        for line in f:
            if line.startswith('#'): continue

            cols       = line.rstrip().split('\t')
            ch         = cols[0]
            curr_class = cols[2]
            key        = (int(cols[3]), int(cols[4]))
            strand     = cols[6]
            
            if curr_class == 'exon' or curr_class == 'mRNA': continue
            else:
                if ch not in gff_dict: gff_dict[ch] = {}
                if strand not in gff_dict[ch]: gff_dict[ch][strand] = OrderedDict() 
                if curr_class not in gff_dict[ch][strand]: gff_dict[ch][strand][curr_class] = OrderedDict() 
                gff_dict[ch][strand][curr_class][key] = ''

def process_matrix():
    """
    Loads gmatrix
    """
    if '.gz' in gmat_fil: f = gzip.open(gmat_fil, 'rb')
    else: f = open(gmat_fil, 'rU')

    for line in f:
        if line.startswith('L'): continue

        locus = line.rstrip().split('\t')[0]
        scaff = locus.split(':')[0]
        pos   = int(locus.split(':')[1])
        if scaff not in snp_dict: snp_dict[scaff] = OrderedDict()
        snp_dict[scaff][pos] = ''
    f.close()

def process_dicts(ch):
    """
    Processes each scaff separately
    """
    dic = OrderedDict()
    try:
        for win,grp in groupby(snp_dict[ch].keys(), lambda x: x // wsize):
            curr_st  = (win * wsize)
            curr_end = ((win + 1) * wsize) - 1
            cds = utr5 = utr3 = intron = inter = 0
            curr_po = list(grp)
            for po in curr_po:
               for strand in gff_dict[ch]:
                   n_intra = sum(1 for i in gff_dict[ch][strand]['gene'].keys() if i[0] <= po <= i[1])
                   if n_intra > 0:
                       n_cds = sum(1 for a in gff_dict[ch][strand]['CDS'].keys() if a[0] <= po <= a[1])
                       cds += n_cds
                       if n_cds == 0:
                           n_5_utr = sum(1 for b in gff_dict[ch][strand]['five_prime_UTR'].keys() if b[0] <= po <= b[1])
                           utr5 += n_5_utr
                           if n_5_utr == 0:
                               n_3_utr = sum(1 for c in gff_dict[ch][strand]['three_prime_UTR'].keys() if c[0] <= po <= c[1])
                               utr3 += n_3_utr
                               if n_3_utr == 0:
                                   intron += 1 
                   else:
                       inter += 1 
            key = (curr_st, curr_end)
            #if ch not in dic: dic[ch] = OrderedDict()
            dic[key] = {'Total': len(curr_po)*2, 'CDS': cds, 
                            "5' UTR": utr5, "3' UTR": utr3, "Intron": intron,
                            "Intergenic": inter}
    except KeyError:
        pass 
    return (ch, dic)

def print_counts(dic_list):
    """
    Prints counts
    """
    ft = ['Intergenic', "5' UTR", "3' UTR", "CDS", "Intron"]
    with open(out_fil, 'wb') as o:
        o.write("CHRM\tST\tEND\tTOTAL\tINTER\t5UTR\t3UTR\tCDS\tINTRON\n")
        for D in sorted(dic_list, key=itemgetter(0)):
            if D[1]:
                for bins in D[1]:
                    #cts = "{Total}\t{Intergenic}\t{5' UTR}\t{3' UTR}\t{CDS}\t{Intron}".format(**D[1][bins])
                    curr_tot = D[1][bins]["Total"]
                    curr_ratios = map("{: 0.4F}".format, [D[1][bins][i] / float(curr_tot) for i in ft])
                    o.write('{0}\t{1[0]}\t{1[1]}\t{2}\t{3}\n'.format(D[0], bins, curr_tot, '\t'.join(curr_ratios)))
    o.close()

if __name__=='__main__':
    start = time.time()
    if len(sys.argv) != 6:
        sys.stderr.write(main.__doc__ + "\n")
        sys.exit()
    gmat_fil = sys.argv[1]
    gff_fil = sys.argv[2]
    out_fil = sys.argv[3]
    n_threads = int(sys.argv[4])
    wsize = int(sys.argv[5])
    main()
    print "Finished; Took:", '{:0.2F}'.format(time.time() - start), "seconds."
