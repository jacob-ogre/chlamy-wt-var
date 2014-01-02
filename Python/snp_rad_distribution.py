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


import sys
import re
import gzip
snp_position_dict = {}
rad_position_dict = {}

pattern = '.{12}GCA.{6}TGC.{12}'
re_obj  = re.compile(pattern)

def main():
    """
    ---------------------------------------------------------------------------
    AUTHOR: Kyle Hernandez
    EMAIL: kmhernan@utexas.edu

    Dataset for histogram of SNP positions along a RAD tag.
    ---------------------------------------------------------------------------
    USAGE: python snp_rad_distribution.py gmatrix.tab reference.fa
    ARGUMENTS:
   	gmatrix.tab  - Tab delimited matrix of snps
 	reference.fa - Fasta file of the Chlamy reference
    """
    process_snps()
    pull_ref_sites()
    print_position_data()

def process_snps():
   try:
       with gzip.open(snpmatrix, 'rb') as f:
            for line in f:
                cols = line.rstrip().split('\t')
                if cols[0] == 'Locus': continue
                else:
                    ch, po = cols[0].split(':')
                    po = int(po)
                    if ch not in snp_position_dict: snp_position_dict[ch] = []
                    snp_position_dict[ch].append(po)
   except:
       with open(snpmatrix, 'rb') as f:
            for line in f:
                cols = line.rstrip().split('\t')
                if cols[0] == 'Locus': continue
                else:
                    ch, po = cols[0].split(':')
                    po = int(po)
                    if ch not in snp_position_dict: snp_position_dict[ch] = []
                    snp_position_dict[ch].append(po)

def pull_ref_sites():
    curr_scaff = ''
    curr_seq   = []
    with open(reference, 'rU') as f:
        for line in f:
            if line[0] == '>' and not curr_scaff:
                curr_scaff = line.rstrip()
            elif line[0] == '>' and curr_scaff:
                joined_seq = ''.join(curr_seq)
                name_scaff = curr_scaff.replace('>', '')
                if name_scaff in snp_position_dict:
                    process_scaff(name_scaff, joined_seq, snp_position_dict[name_scaff])
                curr_scaff = line.rstrip()
                curr_seq = []
            else:
                curr_seq.append(line.rstrip())
        joined_seq =  ''.join(curr_seq)
        name_scaff = curr_scaff.replace('>', '')
        if name_scaff in snp_position_dict:
            process_scaff(name_scaff, joined_seq, snp_position_dict[name_scaff])

def process_scaff(scaff, seq, snp_list):
    cut_sites = re_obj.finditer(seq)
    for site in cut_sites:
        for i in sorted(snp_list):
            if (site.start() + 1) <= i <= site.end():
                curr_snp = i - (site.start() + 1)
                if curr_snp not in rad_position_dict:
                    rad_position_dict[curr_snp] = 0
                rad_position_dict[curr_snp] += 1

def print_position_data():
    print "Position\tCounts"
    for p in sorted(rad_position_dict):
        print str(p) + '\t' + str(rad_position_dict[p])

if __name__=='__main__':
    if len(sys.argv) != 3:
        print main.__doc__
        sys.exit()
    snpmatrix = sys.argv[1]
    reference = sys.argv[2]
    main()
