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
import time
from itertools import groupby, combinations
import gzip
from math import sqrt

def main():
    """
    ---------------------------------------------------------------------------
    AUTHOR: Kyle Hernandez
    EMAIL: kmhernan@utexas.edu

    Return multiple diversity statistics across non-overlapping windows of
    a user-given size.
    ---------------------------------------------------------------------------
    USAGE: python watterson_pi_blockwise.py gmatrix.tab window_size out.tab
    ARGUMENTS:
   	gmatrix.tab - numeric genotype matrix
  	window_size - Size of window in bp (Int)
 	out.tab     - Tab-delimited output file
    """

    snp_table, n = load_gmatrix()
    process_windows(snp_table, n)

def load_gmatrix():
    """
    Loads gmatrix into SNPtable.
    """

    dic = {}
    n = 0
    with gzip.open(gmat, 'rb') as f:
        for line in f:
            if line.startswith('L'):
                n = len(line.rstrip().split('\t')[1::])
                continue

            cols  = line.rstrip().split('\t')
            scf   = cols[0].split(':')[0]
            position = int(cols[0].split(':')[1])
            gtype = cols[1::]

            if gtype.count('N') < 13:
                if scf not in dic: dic[scf] = {}
                dic[scf][position] = gtype
    return dic, n

def process_windows(snp, n):
    """
    Creates and processes windows.
    """
    # Calculate the constant parameters
    # Based on D. Hartl "A Primer of Population Genetics"
    # Theta = S / a1 and Var(Theta) = (theta/k*a1) + (a2 * theta^2 / a1^2)
    # Var(pi)= Var(Theta) = (b1*theta / k) + b2 * theta^2
    # c1 and c2 are used to calculate D
    # k = total number of loci
    a1 = float(sum([1/float(i) for i in xrange(1,n)]))
    a2 = float(sum([1/float(i*i) for i in xrange(1,n)]))
    b1 = (n + 1) / float(3 * (n-1))
    b2 = (2 * ((n * n) + n + 3)) / float(9 * n * (n-1))
    c1 = (b1/a1) - (1 / (a1 * a1))
    c2 = (1 / ((a1 * a1) + a2)) * (b2 - ((n + 2) / (a1 * n)) + (a2 / (a1 * a1)))
    comp = (n * (n - 1)) / 2

    with open(outf, 'wb') as o:
        o.write('Chromosome\tBin\tS\ttheta\tsd_theta\tpi\tsd_pi\tD\tp_typed\tn_loci\n')
        # For each scaffold create windows of w_size bp distance
        for scf in sorted(snp):
            for window, grp in groupby(sorted(snp[scf].keys()), lambda x: x // w_size):
                curr_pos = list(grp)
                n_loci   = len(curr_pos)
                if n_loci > 3:
                    stat_dic = get_stats(curr_pos, snp[scf], n_loci, n, a1, a2,
                                         b1, b2, c1, c2, comp)
                    if stat_dic:
                        fmt_row  = format_stats(scf, n_loci, window, stat_dic)
                        o.write(fmt_row)

def get_stats(positions, dic, n_loci, n_samp, a1, a2, b1, b2, c1, c2, comp):
    """
     Estimates Pi and theta. Based on D. Hartl, A Primer of Population Genetics, 3rd ed.
    """
    # Initialize summary statistics dict
    cts = dict(monomorphic = 0, missingdata = 0, indel = 0, pairwise = 0,
               pi = 0.0, pi_adjust = 0.0, segregate = 0, segregate_adjust = 0.0,
               uninformative = 0.0, segregate_per = 0.0, theta = 0.0,
               theta_adjust = 0.0, tajima = 0.0, tajima_adjust = 0.0,
               denominator = 0.0, denominator_adj = 0.0, p_typed = 0.0,
               sd_theta = 0.0, sd_pi = 0.0)

    # Reduce to polymorphic loci
    align = {}
    typed = 0.0
    for p in sorted(positions):
        typed += 1 - (dic[p].count('N') / float(n_samp))
        if len(set([i for i in dic[p] if i != 'N' and len(i) == 1])) == 1:
            cts["monomorphic"] += 1
        else:
            cts["segregate"] += 1
            for n,i in enumerate(dic[p]):
                if n not in align: align[n] = []
                align[n].append(i)


    if cts["segregate"] > 0:
        # Pairwise data
        for A,B in combinations(sorted(align.keys()), 2):
            seq_A = align[A]
            seq_B = align[B]
            for a,b in zip(seq_A, seq_B):
                if a == 'N' or b == 'N': cts["missingdata"] += 1
                elif len(a) != 1 or len(b) != 1: cts["indel"] += 1
                elif a != b: cts["pairwise"] += 1

        # Check for lots of missing data
        cts["uninformative"] = (cts["missingdata"] + cts["indel"]) / float(n_samp)
        if cts["uninformative"] < int(0.5 * n_loci):
            # Get avg % typed in window and n loci
            cts["p_typed"] = typed / n_loci

            # Get summary statistics and adjustments based on uninformative counts
            # First, estimate S and adjust S by removing uninformative data
            cts["segregate_per"] = cts["segregate"] / float(n_loci)
            cts["segregate_adjust"] = cts["segregate"] / float(n_loci - cts["uninformative"])

            # Estimate theta and adjusted theta all scaled for number of sites
            cts["theta"] = cts["segregate_per"] / a1
            cts["theta_adjust"] = cts["segregate_adjust"] / a1
            cts["sd_theta"] = sqrt((cts["theta"] / (n_loci * a1)) +\
                              ((a2 * (cts["theta"] * cts["theta"])) / \
                              (a1 * a1)))

            # Get pairwise differences to estimate pi
            cts["pi"] = (cts["pairwise"] / float(comp)) / float(n_loci)
            cts["pi_adjust"] = (cts["pairwise"] / float(comp)) /\
                                float(n_loci - cts["uninformative"])
            cts["sd_pi"] = sqrt(((b1 * cts["theta"]) / float(n_loci)) + \
                               (b2 * (cts["theta"] * cts["theta"])))

            # Now estimate Tajima's D
            cts["denominator"] = sqrt((c1 * cts["segregate_per"]) + \
                                     ((c2 * cts["segregate_per"]) * \
                                      (cts["segregate_per"] - (1/float(n_loci)))))
            cts["tajima"] = (cts["pi"] - cts["theta"]) / cts["denominator"]
            cts["tajima_adjust"] = (cts["pi_adjust"] - cts["theta_adjust"]) /\
                                    cts["denominator"]

            return cts

def format_stats(scf, n_loci, window, stat_dic):
    """
    Formats for writing to file
    """
    return '{0}\t{1}\t{2[segregate_per]: 0.5F}\t{2[theta]: 0.5F}\t'.format(scf, window, stat_dic) + \
           '{0[sd_theta]: 0.5F}\t{0[pi]: 0.5F}\t{0[sd_pi]: 0.5F}\t'.format(stat_dic) + \
           '{0[tajima]: 0.5F}\t{0[p_typed]: 0.5F}\t{1}\n'.format(stat_dic, n_loci)

if __name__ == '__main__':
    start = time.time()
    if len(sys.argv) != 4:
        print main.__doc__
        sys.exit()
    gmat   = sys.argv[1]
    w_size = int(sys.argv[2])
    outf   = sys.argv[3]
    main()
    print 'Finished; Took {: 0.3F} seconds.'.format(time.time() - start)
