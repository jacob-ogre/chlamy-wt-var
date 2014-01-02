# SOLiDFilterTrim2bRAD.py
# Filter SOLiD reads per qual scores and trim to 36bp.
#
# CC0, a GPL-compliant copyright.
#
# Written in December 2012 by Jacob Malcom, jacob.malcom@utexas.edu
#
# Although not required in any sense, share the love and pass on attribution
# when using or modifying this code.
#
# To the extent possible under law, the author(s) have dedicated all copyright 
# and related and neighboring rights to this software to the public domain 
# worldwide. This software is distributed without any warranty.
# 
# You should have received a copy of the CC0 Public Domain Dedication along with 
# this software. If not, see <http://creativecommons.org/publicdomain/zero/1.0/>
#

import glob
from itertools import izip
from multiprocessing import Pool
import os
import sys


def main():
    """
    Filter SOLiD reads per qual scores and trim to 36bp.

    USAGE:
        python SOLiDFilterTrim2bRAD.py <tol> <bad> <indir> <outdir>

    ARGS:
        tol, the proportion of bases that are tolerated to be low quality
        bad, the Phred score below which a base is considered low quality
        indir, the path to the directory with raw .csfasta and .qual files
        outdir, the path to which filtered csfasta and qual files are written

    RETURNS:
        files, prefix filt_*, one for the .csfasta and one for the .qual inputs

    COMMENTS:
        tol and bad should be set at a level with which the researcher is 
        comfortable; I prefer tol=0.25 and bad=16 (~p=0.05) at this time, but
        more stringent requirements may reduce false-positive SNPs.

        NOTE that the number of processors for multiprocessing is hard-coded to
        12, which is suitable for Lonestar and my MacPro, but may need to be
        changed if used on a different system.

        UPDATE: I was, originally, translating the first cs base of the read,
        then trimming the head and returning to position 38.  It appears, 
        however, that SHRiMP really does assume the adaptor base is present, and
        I need to actually keep the first (letterspace) base and trim at 37.
    """
    os.chdir(indir)
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    seq_files = glob.glob("*.csfasta")

    pool = Pool(processes=12)
    pool.map(trim_filter, seq_files)

def trim_filter(fil):
    """
    Process sequence and qual files in parallel to filter.

    ARGS:
        fil, a .csfasta file 
    RETURNS:
        two files, one each filtered and trimmed .csfasta and .qual
    """
    sample = fil.split('.')[0]
    qua_fil = sample + "_QV.qual"
    suffix = '_' + sample[-2:]

    # open the output (filtered) files:
    fil_seq = open(outdir + "/filt_" + fil, "wb")
    fil_qua = open(outdir + "/filt_" + qua_fil, "wb")

    # open the input files:
    seq = open(fil)
    qua = open(qua_fil)

    cur_read, cur_seq, cur_qua = '', '', ''
    for line1, line2 in izip(seq, qua):
        if line1[0] == '>' and cur_read == '':
            if line1.rstrip() != line2.rstrip():
                print "Not same order!\n", line1, line2
            cur_read = line1.rstrip() + suffix + '\n'
        elif line1[0] == '>' and cur_read != '':
            if test_good(cur_seq, cur_qua):
                fil_seq.write(cur_read)
                fil_seq.write(cur_seq)
                fil_qua.write(cur_read)
                qua_line = ' '.join(cur_qua) + '\n'
                fil_qua.write(qua_line)
            if line1.rstrip() != line2.rstrip():
                print "Not same order!\n", line1, line2
            cur_read = line1.rstrip() + suffix + '\n'
        else:
            cur_seq = trim(line1)
            cur_qua = trim(line2)
    fil_seq.close()
    fil_qua.close()
    seq.close()
    qua.close()

def test_good(seq, qua):
    """Basic tests to determine if a read is high-quality."""
    # quickly reject samples with Ns (.)
    if "." in seq:
        return False

    # tests if the proportion of qual scores is below the tolerance:
    quals = [int(x) for x in qua]
    if (sum([x < bad for x in quals]) / 36.0) < tol:
        return True
    else:
        return False

def trim(line):
    """
    Trims input line to len==37.

    ARGS:
        line, which may be csfasta sequence or qual string
    RETURNS:
        if line is csfasta, returns a trimmed string
        if line is qual, returns a split list to be joined later for writing
    """
    if len(line) < 55:
        return line[0:37] + '\n'
    else:
        return line.split(' ')[0:36]

if __name__ == '__main__':
    if len(sys.argv) < 5:
        print main.__doc__
        sys.exit()

    tol = float(sys.argv[1])
    bad = int(sys.argv[2])
    indir = sys.argv[3]
    outdir = sys.argv[4]
    main()