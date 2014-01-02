# VCFHighQualityVariants.py
# Extract high-quality variants from .vcf and save to new .vcf.
#
# CC0, a GPL-compliant copyright.
#
# Written in January 2013 by Jacob Malcom, jacob.malcom@utexas.edu
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

import os
from scipy import stats 
import sys

def main():
    """
    Filter .vcf files to remove all but very high-quality SNPs for BQRC.

    USAGE: 
        python VCFHighQualityVariants.py <infile> <outdir>
    
    ARGS:
        infile, the name of the .vcf file to be filtered
        outdir, the name of the directory to which the qual files are written
    
    RETURNS:
        A filtered file with name *_HQ_SNPs.vcf, using the base (line number) of
        infile as the prefix.
    """
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    MQs = extract_quals(infile)
    cutoff = stats.scoreatpercentile(MQs, 90)
    filter_file(cutoff)

def extract_quals(f):
    """Extract MAPQs for variant loci in file f; return list"""
    MQ = []
    for line in open(f):
        if not line.startswith('#'):
            split = line.rstrip().split('\t')
            if split[4] != '.':
                MQ.append(float(split[5]))
    return MQ

def filter_file(cutoff):
    """Writes variants with GQ > 90 and MAPQ > 90th percentile (from tol)."""
    base = infile.split('.')[0].split('/')[-1].split('_')[0]
    print base
    outfil = outdir + '/' + base + "_HQ_SNPs.vcf"
    with open(outfil, 'wb') as out:
        for line in open(infile):
            if line.startswith("#"):
                out.write(line)
            else:
                split = line.rstrip().split('\t')
                if split[4] != '.':
                    geno_dat = split[9].split(':')
                    if float(geno_dat[3]) >= 90 and float(split[5]) > cutoff:
                        out.write(line)


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print main.__doc__
        sys.exit()

    infile = sys.argv[1]
    outdir = sys.argv[2]
    main()

