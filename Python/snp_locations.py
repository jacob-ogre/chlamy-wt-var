#!/usr/bin/env python
# Calculate the distribution of polymorphic RAD loci across site classes.

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

gff_dict  = {}
snp_dict  = {}
gene_dict = {}
def main():
    """
    ---------------------------------------------------------------------------
    AUTHOR: Kyle Hernandez
    EMAIL: kmhernan@utexas.edu

    Calculate the distribution of polymorphic RAD loci across site classes.
    ---------------------------------------------------------------------------

    USAGE: python snp_locations.py gmatrix.tab file.gff out.tab n_threads

    ARGUMENTS:
    	gmatrix.tab - Tab-delimited genotype matrix file of variant sites
        file.gff    - GFF file
        out.tab     - Output file of counts
        n_threads   - The number of threads to run
    """

    # Load the GFF and SNP positions into dictionaries
    load_gff()
    intergenic = process_matrix()
    
    # Map:
    # Create a pool of n_threads workers and use them to process
    # scaffolds separately
    ch_vals = sorted(gff_dict.keys())
    sys.stdout.write("Counting features...\n")
    pool    = Pool(processes=n_threads)
    ct_list = pool.map(process_dicts, ch_vals)

    # Reduce:
    # Process the list of dicts
    print_counts(intergenic, ct_list)

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
            if curr_class == 'gene':
                if ch not in gene_dict: gene_dict[ch] = {}
                if strand not in gene_dict[ch]: gene_dict[ch][strand] = {}
                gene_dict[ch][strand][key] = ''

            else:
                if ch not in gff_dict: gff_dict[ch] = {}
                if strand not in gff_dict[ch]: gff_dict[ch][strand] = {} 
                if curr_class not in gff_dict[ch][strand]: gff_dict[ch][strand][curr_class] = {} 
                gff_dict[ch][strand][curr_class][key] = ''

def process_matrix():
    sys.stdout.write("Loading Genotype Matrix dictionary...\n")
    intergenic = 0
    tot = 0
    with gzip.open(gmat_fil, 'rb') as f:
        for line in f:
            if line.startswith('Locus'): continue
            cols = line.rstrip().split('\t')
            ch   = cols[0].split(':')[0]
            po   = int(cols[0].split(':')[1])
            tot += 1
            try:
                curr_d = gene_dict[ch]
                for strand in gene_dict[ch]:
                    d = []
                    while not d:
                        d = [k for k in gene_dict[ch][strand].keys() if k[0]<= po <= k[1]]
                        break
                    if d:
                        if ch not in snp_dict: snp_dict[ch] = {}
                        snp_dict[ch][po] = ''
                    else:
                        intergenic += 1
                    
            except KeyError:
                continue
    print "Total SNPs:", tot
    return intergenic

def process_dicts(ch):
    """
    Process each scaffold separately, return dict of counts.
    """

    try:
        dic = {}
        for po in sorted(snp_dict[ch].keys()):
            for strand in sorted(gff_dict[ch]):
                curr_ch = gff_dict[ch][strand].keys()
                n_cds = len([a for a in gff_dict[ch][strand]['CDS'].keys() if a[0] <= po <= a[1]])
                if n_cds > 0:
                    if 'CDS' not in dic: dic['CDS'] = 0
                    dic['CDS'] += n_cds
                else:
                    n_5_utr = len([b for b in gff_dict[ch][strand]['five_prime_UTR'].keys() if b[0] <= po <= b[1]])
                    if n_5_utr > 0:
                        if "5' UTR" not in dic: dic["5' UTR"] = 0
                        dic["5' UTR"] += n_5_utr
                    else:
                        n_3_utr = len([c for c in gff_dict[ch][strand]['three_prime_UTR'].keys() if c[0] <= po <= c[1]])
                        if n_3_utr > 0:
                            if "3' UTR" not in dic: dic["3' UTR"] = 0
                            dic["3' UTR"] += n_3_utr
                        else:
                            if "Intron" not in dic: dic["Intron"] = 0
                            dic["Intron"] += 1 

    except: pass 
    return dic

def print_counts(intergenic, dic_list):
    """
    Reduce and print the list of dicts.
    """
    features = sorted(dic_list[0].keys())
    
    with open(out_fil, 'wb') as o:
        o.write('Class\tCounts\tSum\tRatio\n')

        # List comprehension to sum all counts
        total_sum = sum(d[f] for d in dic_list \
                        for f in features if d and \
                        f in d) + intergenic
        o.write("Intergenic\t{0}\t{1}\t{2: 0.3F}\n".format(intergenic, total_sum, intergenic/float(total_sum)))
        
        for f in features:
            # List comprehension to sum counts for feature f
            curr_sum = sum(d[f] for d in dic_list \
                           if d and f in d)
            curr_ratio = curr_sum/float(total_sum)

            o.write("{0}\t{1}\t{2}\t{3: 0.3F}\n".format(f, curr_sum, total_sum, curr_ratio))

if __name__=='__main__':
    start = time.time()
    if len(sys.argv) != 5:
        sys.stderr.write(main.__doc__ + "\n")
        sys.exit()
    gmat_fil = sys.argv[1]
    gff_fil = sys.argv[2]
    out_fil = sys.argv[3]
    n_threads = int(sys.argv[4])
    main()
    print "Finished; Took:", '{:0.2F}'.format(time.time() - start), "seconds."
