#!/usr/bin/env python
# Smooth over moving windows.
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
from itertools import groupby, izip

def main():
    """
    ---------------------------------------------------------------------------
    AUTHOR: Kyle Hernandez
    EMAIL: kmhernan@utexas.edu

    Averages over nwin windows
    ---------------------------------------------------------------------------
    USAGE: python smooth_windows.py in.tab wsize nwin out.tab
    ARGUMENTS:
	in.tab  - Input tab-delimited file from watterson_pi_window.py
        wsize   - Size of windows
   	nwin    - Number of adjacent windows to smooth over
 	out.tab - Output tab-delimited file
    """
    stat_dic, head = load_stats()
    smooth_windows(stat_dic, head)

def load_stats():
    """
    Loads in.tab into dic
    """
    dic = {}
    head = ['Chromosome', 'Start', 'Stop']
    with open(dat, 'rU') as f:
        for line in f:
            if line.startswith('Chromosome'): 
                [head.append(j) for j in line.rstrip().split('\t')[2::]]
                continue
            if line.startswith('chromosome'):
                cols = [float(i) if n > 1 else int(i) if n == 1 or n == 9 
                        else i for n,i in enumerate(line.rstrip().split('\t'))]
                ch = cols[0]
                bn = ((cols[1] + 1) * wsize) - 1
                if ch not in dic: dic[ch] = {}
                dic[ch][bn] = cols[2::]
    return dic, head

def smooth_windows(stat, head):
    """
    Takes averages over nwin windows
    """
    adjust = wsize * nwin
    with open(out, 'wb') as o:
        o.write('\t'.join(head) + '\n')
        for ch in sorted(stat):
            for b,grp in groupby(sorted(stat[ch].keys()), lambda x: x // adjust):
                upper = ((b + 1) * adjust) - 1 
                lower = b * adjust
                curr_bins = list(grp)
                curr_dat = [stat[ch][i] for i in curr_bins]
                curr_avg = [sum(x)/float(len(x)) for x in izip(*curr_dat)]
                o.write(get_row(ch, lower, upper, curr_avg))

def get_row(ch, lower, upper, data):
    """
    Formats means 
    """
    return '{0}\t{1}\t{2}\t{3[0]: 0.5F}\t{3[1]: 0.5F}\t{3[2]: 0.5F}\t'.format(ch, lower, upper, data) + \
           '{0[3]: 0.5F}\t{0[4]: 0.5F}\t{0[5]: 0.5F}\t{0[6]: 0.5F}\t{0[7]: 0.2F}\n'.format(data)

def grouper(n, iterable, fillvalue=None):
    """
    Collect data into fixed-length chunks or blocks.
    Taken from python cookbook.
    """
    args = [iter(iterable)] * n
    return izip_longest(fillvalue=fillvalue, *args)
 
if __name__ == '__main__':
    start = time.time()
    if len(sys.argv) != 5:
        print main.__doc__
        sys.exit()
    dat = sys.argv[1]
    wsize = int(sys.argv[2])
    nwin = int(sys.argv[3])
    out = sys.argv[4]
    main()
    print 'Finished; Took: {0: 0.4F} seconds.'.format(time.time() - start)
