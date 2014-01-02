# VCFGenotypeQualityDistribution.py
# Extract GQ for ALT alleles from a .vcf file to look at histogram.
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

import glob
import os
import sys

def main():
    """
    Filter .vcf files to remove all but very high-quality SNPs for BQRC.

    USAGE: 
        python FilterSAM.py <tol> <outdir>
    
    ARGS:
        outdir, the name of the directory to which the qual files are written
    
    RETURNS:
        .csv file of with cols Type, MAPQ, GQ
    """
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    for fil in glob.glob("*GATK.vcf"):
        qual_dict = extract_quals(fil)
        print_res(fil, qual_dict)

def extract_quals(f):
    dic = {"gqi": [], 
           "gqs": [], 
           "MQ": []}
    for line in open(f):
        if not line.startswith('#'):
            split = line.rstrip().split('\t')
            if split[4] != '.':
                dic["MQ"].append(split[5])
                geno_dat = split[9].split(':')
                if len(split[4]) == 1:
                    dic["gqs"].append(geno_dat[3])
                else:
                    dic["gqi"].append(geno_dat[3])
    return dic

def print_res(f, dic):
    new = outdir + '/' + f.split('_')[0] + "_varQuals.csv"
    GQs = dic["gqs"] + dic["gqi"]
    gqs_ind = ("SNP " * len(dic["gqs"])).split(' ')[:-1]
    gqi_ind = ("indel " * len(dic["gqi"])).split(' ')[:-1]
    types = gqs_ind + gqi_ind
    with open(new, "wb") as out:
        header = "Type" + ',' + "MAPQ" + ',' + "GQ" + '\n'
        out.write(header)
        for i in range(len(GQs)):
            line = types[i] + ',' + dic["MQ"][i] + ',' + GQs[i] + '\n'
            out.write(line)

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print main.__doc__
        sys.exit()

    outdir = sys.argv[1]
    main()

