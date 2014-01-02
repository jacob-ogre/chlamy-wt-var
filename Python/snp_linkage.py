# Calculate LD between SNPs.
# Copyright (C) 2014 K. Hernandez

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


import sys
import scipy.stats
import time
from itertools import combinations, groupby
import numpy as np
from operator import itemgetter
import os

def main():
    """
    ---------------------------------------------------------------------------
    AUTHOR: Kyle Hernandez
    EMAIL: kmhernan@utexas.edu

    Estimates correlations between all pairwise positions on a scaffold.
    ---------------------------------------------------------------------------
  
    USAGE: python snp_linkage.py gmatrix.tab chl_ref.fa est_type max_dist window_cutoff

    ARGUMENTS:
    	gmatrix.tab   - Genotype matrix file
        ch_ref.fa     - Chlamy reference fasta file
        est_type      - Estimate across all scaffolds or within each [all/within] 
        max_dist      - Maximum distance to consider [Int]
        window_cutoff - Value separating moving vs. static windows. [Int]

    RETURNS:
  	Automatically outputs to file in ~/Projects/Chlamy/RAD/results/WT_RAD/LD_data with
        the est_type, max_dist, and window_cutoff values in filename
    """
    snp_data = load_snps()
    process_snps(snp_data)

def load_snps():
    """
    Loads the genotype matrix SNP data into a dictionary
    """
    dic = {}
    multi_allelic = 0
    low_coverage  = 0
    with open(gmatrix, 'rU') as f:
        for line in f:
            if line.startswith('Locus'): continue

            cols = line.rstrip().split('\t')
            chrm = cols[0].split(':')[0]
            pos  = int(cols[0].split(':')[1])
            calls = cols[1::]
            if len(set([i for i in calls if i != 'N'])) > 2:
                multi_allelic += 1
            elif len(calls) - calls.count('N') < 3:
                low_coverage += 1
            else:
                if chrm not in dic: dic[chrm] = {}
                dic[chrm][pos] = calls
    print "# Loci > 2 Alleles:", multi_allelic
    print "# Loci with < 3 typed individuals:", low_coverage
    return dic

def process_snps(snp_data):
    """
    Loads the chlamy reference and checks whether the current scaffold is in the dictionary
    of snp_data. Passes the scaffold to the estimate_LD function.
    """
    curr_scaff = ''
    curr_seq   = []
    o_fil      = odir + os.sep + 'CR_' + est_type + '_MD' +\
                 str(max_dist) + '_WC' + str(window_cutoff) + '.tab'
    o = open(o_fil, 'wb')

    if est_type == 'within':
        o.write('CHRM\tDIST\tRSQ\n')
    elif est_type == 'all':
        o.write('DIST\tRSQ\n')
        all_D_dict = {} 

    print 'Processing reference...'
    with open(ch_ref, 'rU') as f:
        for line in f:

            if line.startswith('>') and not curr_scaff:
                curr_scaff = line.rstrip()

            elif line.startswith('>') and curr_scaff:
                joined_seq = ''.join(curr_seq)
                name_scaff = curr_scaff.replace('>','')

                if est_type == 'within':
                    if name_scaff in snp_data:
                        curr_positions = snp_data[name_scaff]
                        within_LD(name_scaff, curr_positions, joined_seq, o)
                elif est_type == 'all':
                    if name_scaff in snp_data:
                        curr_positions = snp_data[name_scaff]
                        all_D_dict = all_LD(curr_positions, joined_seq, all_D_dict)
                    
                curr_scaff = line.rstrip()
                curr_seq   = []

            else:
                curr_seq.append(line.rstrip())

        joined_seq = ''.join(curr_seq)
        name_scaff = curr_scaff.replace('>', '')
        if est_type == 'within':
            if name_scaff in snp_data:
                curr_positions = snp_data[name_scaff]
                within_LD(name_scaff, curr_positions, joined_seq, o)
        elif est_type == 'all':
            if name_scaff in snp_data:
                curr_positions = snp_data[name_scaff]
                all_D_dict = all_LD(curr_positions, joined_seq, all_D_dict)

    if est_type == 'all':
        process_all_type(all_D_dict, o)

    o.close()

def within_LD(scaffold, positions, sequence, o):
    """
    Estimate LD within a scaffold over a moving-average window based on the user
    input values: max_dist and window_cutoff 
    """
    D_dict = {}
    for i,j in combinations(sorted(positions.keys()), 2):
        distance = abs(i-j)
        if distance <= max_dist:
            if distance not in D_dict: D_dict[distance] = []
            D_dict[distance].append(get_rho(sequence[i-1], sequence[j-1],\
                                       positions[i], positions[j]))
    for d in sorted(D_dict):
        if d <= window_cutoff:
            curr_range = range(d, (d+1) * 2)
            curr_values = []
            for i in [k for k in curr_range if k in D_dict]:
                [curr_values.append(l) for l in D_dict[i]]
            if len(curr_values) > 3:
                curr_mean = np.mean(np.array(curr_values))
                o.write(scaffold + '\t' + str(d) + '\t' + '{:.4F}'.format(curr_mean) + '\n')

        else:
            curr_range = range(d, (d + (window_cutoff * 2) + 2))
            curr_values = []
            for i in [k for k in curr_range if k in D_dict]:
                [curr_values.append(l) for l in D_dict[i]]
            if len(curr_values) > 3:
                curr_mean = np.mean(np.array(curr_values))
                o.write(scaffold + '\t' + str(d) + '\t' + '{:.4F}'.format(curr_mean) + '\n')
                
    print 'Finished with', scaffold

def all_LD(positions, sequence, dic):
    for i,j in combinations(sorted(positions.keys()), 2):
        distance = abs(i-j)
        if distance <= max_dist:
            if distance not in dic: dic[distance] = []
            dic[distance].append(get_rho(sequence[i-1], sequence[j-1],\
                                       positions[i], positions[j]))

    return dic

def get_rho(ref_A, ref_B, gt_A, gt_B):
    A = np.array([0 if k == ref_A else 3 if k == 'N' \
                  else 1 for k in gt_A],\
                  dtype=np.float)
    B = np.array([0 if l == ref_B else 3 if l == 'N' \
                  else 1 for l in gt_B],\
                  dtype=np.float)
    A[A==3.]=np.nan
    B[B==3.]=np.nan
    rho, pval = scipy.stats.spearmanr(A, B)
    if not np.isnan(rho):
        return (rho ** 2)

def process_all_type(D_dict, o):
    """
    Calculates moving window rhosquared average across all chromosomes.
    Prints to file.
    """
    print 'Writing to file...'
    for d in sorted(D_dict):
        if d <= window_cutoff:
            curr_range = range(d, (d+1) * 2)
            curr_values = []
            for i in [k for k in curr_range if k in D_dict]:
                [curr_values.append(l) for l in D_dict[i]]
            if len(curr_values) > 3:
                curr_mean = np.mean(np.array(curr_values))
                o.write(str(d) + '\t' + '{:.4F}'.format(curr_mean) + '\n')

        else:
            curr_range = range(d, (d + (window_cutoff * 2) + 2))
            curr_values = []
            for i in [k for k in curr_range if k in D_dict]:
                [curr_values.append(l) for l in D_dict[i]]
            if len(curr_values) > 3:
                curr_mean = np.mean(np.array(curr_values))
                o.write(str(d) + '\t' + '{:.4F}'.format(curr_mean) + '\n')
 
if __name__=='__main__':
    start = time.time()
    VALID_TYPES = ['within', 'all']

    if len(sys.argv) != 7:
        print main.__doc__
        sys.exit()
    elif sys.argv[3] not in VALID_TYPES:
        print "ERROR! Incorrect option for est_type!!!"
        print main.__doc__
        sys.exit()

    gmatrix = sys.argv[1]
    ch_ref  = sys.argv[2]
    est_type = sys.argv[3]
    max_dist = int(sys.argv[4])
    window_cutoff = int(sys.argv[5])
    odir = os.path.abspath(sys.argv[6])
    main()
    print "Finished; Took:", '{:5.5F}'.format(time.time() - start), "seconds."
