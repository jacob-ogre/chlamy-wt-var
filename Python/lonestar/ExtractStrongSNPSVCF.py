# ExtractStrongSNPsVCF.py
# Filter .vcf files to remove all but very high-quality SNPs for BQRC.
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
import os
import sys

def main():
    """
    Filter .vcf files to remove all but very high-quality SNPs for BQRC.

    USAGE: 
        python FilterSAM.py <tol> <cut_fil> <outdir>
    
    ARGS:
        tol, MAPQ lower bound, MAPQ > tol (not >=)
        cut_fil, a .csv of 90th percentile MAPQs for each line
        outdir, the name of the directory to which the filtered files are 
            written
    
    RETURNS:
        Filtered .vcf containing only high-quality variants.
    
    COMMENTS:
        Globs for files on the pattern *GATK.vcf in the pwd.
    """
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    MAPQ_dic = get_mapq_dat(cut_fil)

    for fil in glob.glob("*GATK.vcf"):
        filter_file(fil, MAPQ_dic)


def get_mapq_dat(fil):
    """Extract 90th percentile cutoff from .csv fil and return dic."""
    dic = {}
    for line in open(fil):
        if not line.startswith('\"\"'):
            new = line.replace('"', '').rstrip().split(',')
            dic[new[1]] = float(new[2])
    return dic


def filter_file(fil, dic):
    """Writes variants with GQ > tol and MAPQ > 90th percentile from dic."""
    base = fil.split('.')[0].split('_')[0]
    cutoff = dic[base]
    outfil = outdir + '/' + base + "_HQ_SNPs.vcf"
    n_retained = 0
    with open(outfil, 'wb') as out:
        for line in open(fil):
            if line.startswith("#"):
                out.write(line)
            else:
                split = line.rstrip().split('\t')
                if split[4] != '.':
                    geno_dat = split[9].split(':')
                    if float(geno_dat[3]) >= tol and float(split[5]) > cutoff:
                        out.write(line)
                        n_retained += 1
    print fil, n_retained


if __name__ == '__main__':
    if len(sys.argv) < 4:
        print main.__doc__
        sys.exit()

    tol = int(sys.argv[1])
    cut_fil = sys.argv[2]
    outdir = sys.argv[3]
    main()