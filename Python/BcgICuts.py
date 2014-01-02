# BcgICuts.py
# Finds, tabulates, and writes all BcgI cut sites in a genome file.
# Copyright (C) 2013 Jacob Malcom, jmalcom@uconn.edu

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


import re
import os
import sys
import time
# import pprint
from operator import ne
from itertools import imap, combinations
import numpy as np
# from scipy.spatial import distance

bcg_dic = {}
uniq_seqs = {}
bcgi_uniq = {}

def main():
    """
    Finds, tabulates, and writes all AlfI cut sites in a genome file.

    USAGE:
        python BcgICuts.py <infil> <outfil>
    
    ARGS:
        <infil>, a fasta-formatted sequence file

    RETURNS:
        <outfil>, a fasta file with one entry per cut site
        A summary of the process (# cut sites, etc.)

    COMMENTS:
        Uses the re module to compile and search for the regexes 
            '[ACGTN]{12}CGA[ACGTN]{6}TGC[ACGTN]{12}' and
            '[ACGTN]{12}GCA[ACGTN]{6}TCG[ACGTN]{12}'
        which corresponds to the BcgI seq. After finding all sites, winnows the 
        set to those with a unique recognition sequence (mapping cannot 
        differentiate sites if identical sequence).

    NOTE:
        This script includes commented-out code for calculating Hamming distance
        between pairs of cutsites. This process is quite time-consuming, so the
        distances are simply sampled 1 of every 10000 combinations, which still
        gives a representative look at the distribution of distances between
        sequences, and thus an idea of how unique sequences are around AlfI
        sites.

        Also note that I should have written this not to use bcg_dic, uniq_seqs,
        and bcgi_uniq as global variables...
    """
    bcg_re1 = compile_re1()
    bcg_re2 = compile_re2()
    process_file(infil, bcg_re1, bcg_re2)
    winnow_to_unique()
    write_dict(bcgi_uniq, outfil)
    # dist = check_distances()
    # write_diff(dist, "Chlamy_HammingD.out")
    print_summary()

def process_file(fil, re1, re2):
    """
    Join lines in FASTA to single line, search with recomps to find cut sites.

    Scaffolds in the reference need to be searched in-whole, i.e., in case there
    are recognition sites that span eols. ''.join(list) is much faster and 
    memory efficient than concatenating the strings together. Note that a final 
    call to join the seq, finditer, and add_to_dict is required to process the
    final scaffold in the file.
    """
    current_name = ''
    current_seq = []
    with open(fil) as f:
        for line in f:
            if line[0] == '>' and current_name == '':
                current_name = line.rstrip()
            elif line[0] == '>' and current_name != '':
                joined_seq = ''.join(current_seq)
                plus_sites = re1.finditer(joined_seq)
                neg_sites = re2.finditer(joined_seq)
                add_to_dict(current_name, plus_sites)
                add_to_dict(current_name, neg_sites)
                current_name = line.rstrip()
                current_seq = []
            else:
                current_seq.append(line.rstrip())
        joined_seq = ''.join(current_seq)
        plus_sites = re1.finditer(joined_seq)
        neg_sites = re2.finditer(joined_seq)
        add_to_dict(current_name, plus_sites)
        add_to_dict(current_name, neg_sites)

def add_to_dict(scaf, sites):
    """Parses re MatchObj iterator to add to dict."""
    for site in sites:
        key = ''.join([scaf, ':', str(site.start()), '-', str(site.end())])
        bcg_dic[key] = site.group()     # 'group' is the sequence

def winnow_to_unique():
    """
    Determines if the sequence at a locus is unique.
    
    First goes through alf_dic and counts the number of times each unique
    sequence is found.  Next adds loci to alfi_uniq iff count == 1.
    """
    for key in bcg_dic:
        if bcg_dic[key] not in uniq_seqs:
            uniq_seqs[bcg_dic[key]] = 1
        else:
            uniq_seqs[bcg_dic[key]] += 1
    for key in bcg_dic:
        if uniq_seqs[bcg_dic[key]] == 1:
            bcgi_uniq[key] = bcg_dic[key]

def check_distances():
    """Use itertools combinations for unique pairs; return list of distances."""
    distances = []
    seqs = alfi_uniq.values()
    for i, j in combinations(seqs, 2):
        distances.append(hamming(i, j))
    return distances

def hamming(i, j):
    """Return the number of differences between seqs i, j."""
    return sum(imap(ne, i, j))

def print_summary():
    print "Total # sites:", len(bcg_dic)
    print "# unique sites", len(bcgi_uniq)

    # dic = {}
    # dic.fromkeys(distances, 0)
    # for k in distances:
    #     if k not in dic:
    #         dic[k] = 1
    #     else:
    #         dic[k] += 1
    # print "Hamming distances:"
    # pprint.pprint(dic)


def compile_re1():
    """Compiles the BcgI regex for + strand."""
    pattern = '.{12}CGA.{6}TGC.{12}'
    re_obj = re.compile(pattern)
    return re_obj

def compile_re2():
    """Compiles the BcgI regex for - strand."""
    pattern = '.{12}GCA.{6}TCG.{12}'
    re_obj = re.compile(pattern)
    return re_obj   

# def write_diff(distances, fil):
#   with open(fil, 'wb') as out:
#       for i in distances:
#           line = str(distances[i]) + '\n'
#           out.write(line)

def write_dict(dic, fil):
    with open(fil, 'wb') as f:
        for k in dic:
            line = ''.join([k, '\n', dic[k], '\n'])
            f.write(line)


if __name__ == '__main__':
    if len(sys.argv) < 3:
        print main.__doc__
        sys.exit()

    infil = sys.argv[1]
    outfil = sys.argv[2]

    start = time.time()
    main()
    print "Time required:", time.time() - start, "seconds"
