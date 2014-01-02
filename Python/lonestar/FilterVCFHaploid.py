# FilterVCF.py
# Filter .vcf files to remove low-quality or low-coverage typed loci.
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

# import glob
import os
import sys

def main():
    """
    Filter .vcf files to remove low-quality typed loci.

    USAGE: 
        python FilterSAM.py <tol> <outdir>
    
    ARGS:
        tol is the minimum Phred score of a genotype call
        outdir, the name of the directory to which the filtered files are 
            written
    
    RETURNS:
        Filtered files of name *_filtered.vcf 
    
    COMMENTS:
        Globs for files on the pattern *.vcf in the pwd.
    """
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    # for fil in glob.glob("*SAM.vcf"):
    #     summary = filter_file(fil)
    #     print_results(fil, summary)
    summary = filter_file(infil)
    print_results(infil, summary)

def filter_file(fil):
    res = {"n_typed": 0,
           "n_poly": 0,
           "n_indel": 0,
           "n_SNP": 0,
           "n_invar": 0,
           "n_het": 0}
    outfil = outdir + '/' + fil.split('.')[0] + "_filtered.vcf"
    with open(outfil, 'wb') as out:
        for line in open(fil):
            if line.startswith("#"):
                out.write(line)
            else:
                split = line.rstrip().split('\t')
                if ':' not in split[9]:
                    # set at p > 0.01 for keeping a ref allele:
                    if float(split[9]) < 20 and float(split[5]) > tol:
                        res["n_typed"] += 1
                        res["n_invar"] += 1
                        out.write(line)
                else:
                    geno_dat = split[9].split(':')
                    if float(geno_dat[2]) > tol:
                        likelihoods = [int(x) for x in geno_dat[1].split(',')]
                        try:
                            # these indices correspond to the diagonal of the
                            # matrix of possible diploid genotypes (up to four
                            # alleles), which indicates homs:
                            if likelihoods.index(min(likelihoods)) == 0 or \
                                likelihoods.index(min(likelihoods)) == 2 or \
                                likelihoods.index(min(likelihoods)) == 5 or \
                                likelihoods.index(min(likelihoods)) == 9:
                                res["n_poly"] += 1
                                out.write(line)
                                if len(split[3]) == 1 and len(split[4]):
                                    res["n_SNP"] += 1
                                else:
                                    res["n_indel"] += 1
                            else:
                                # here it will get a little tricky to extract 
                                # the genotypes, but casting the likelihoods to
                                # an array should be pretty easy.
                                res["n_het"] += 1
                        except:
                            print line
    return res

def print_results(fil, dic):
    """Print the summary data to stdout."""

    # Create a few more dict keys & values
    dic["prop_het"] = dic["n_het"] / float(dic["n_typed"])
    dic["p_ind_all"] = dic["n_indel"] / float(dic["n_typed"])
    dic["p_ind_poly"] = dic["n_indel"] / float(dic["n_poly"])
    dic["p_SNP_all"] = dic["n_SNP"] / float(dic["n_typed"])
    dic["p_SNP_poly"] = dic["n_SNP"] / float(dic["n_poly"])

    results = """
# genotyped loci:\t{n_typed}
# polymorphic:\t{n_poly}
# indels:\t{n_indel}
# SNPs:\t{n_SNP}
# hets:\t{n_het}
prop. hets:\t{prop_het}
prop. indel (all):\t{p_ind_all}
prop. indel (poly):\t{p_ind_poly}
prop. SNP (all):\t{p_SNP_all}
prop. SNP (poly):\t{p_SNP_poly}
"""

    print "\n", fil
    print "Given tol={0}".format(tol)
    print results.format(**dic)

if __name__ == '__main__':
    if len(sys.argv) < 3:
        print main.__doc__
        sys.exit()

    tol = int(sys.argv[1])
    infil = sys.argv[2]
    outdir = sys.argv[3]
    main()
