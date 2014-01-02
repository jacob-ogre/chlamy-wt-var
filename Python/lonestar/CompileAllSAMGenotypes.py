# CompileAllSAMGenotypes.py
# Compile genotypes from multiple .vcf files into line x genotype .tab files.
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
    Compile genotypes from multiple .vcf files into line x genotype .tab files.

    USAGE:
        python CompileAllGenotypes.py <allfil> <polyfil> <propcut> <propfil> 
                                      <alltyped>

    ARGS:
        allfil, the name of the file of all typed loci for all lines
        polyfil, the name of the file of all polymorphic loci for all lines
        propcut, the proportion of lines that must be typed at a locus for the 
            locus to be added to propfil
        propfil, the name of the file of polymorphic loci for all lines where 
            the propotion of typed lines is > propcut
        alltyped, the name of the file of loci with all lines typed

    RETURNS:
        four files, allfil, polyfil, propfil, and alltyped; tab-delimited

    COMMENTS:
        Uses glob on the pattern *filtered.vcf to pull in .vcf files; note that
        the filtered.vcf files have already had low-quality genotype calls 
        removed.  Although each of the write_* functions could be cleaned up to 
        shorten the script, it's still pretty fast (~5 min currently).  
    """
    os.chdir(os.getcwd())

    all_loci = {}
    fil_list = glob.glob("*.vcf")
    for fil in fil_list:
        all_loci = get_all_loci(fil, all_loci)

    for fil in fil_list:
        all_loci = fill_dict(fil, all_loci)

    write_all(all_loci, fil_list)
    write_poly(all_loci, fil_list)
    write_prop(all_loci, fil_list)
    write_alltyped(all_loci, fil_list)

def get_all_loci(f, d):
    """Read loci from file f and add to dict d."""
    for line in open(f):
        if not line.startswith('#'):
            split = line.rstrip().split('\t')
            if ':' in split[9]:
                GQ = float(split[9].split(':')[-1])
            else:
                # note that split[9] is PL, not GQ, but QUAL works here:
                GQ = float(split[5])
            if float(split[5]) > 30 and GQ > 30:
                chrom = split[0]
                locus = int(split[1])
                if chrom in d:
                    d[chrom][locus] = []
                else:
                    d[chrom] = {locus: []}
    return d

def fill_dict(f, d):
    """Fill dict d with genotypes from file f if typed; else N."""
    dic = {}
    for line in open(f):
        if not line.startswith('#'):
            split = line.split('\t')
            if ':' in split[9]:
                GQ = float(split[9].split(':')[-1])
            else:
                # note that split[9] is PL, not GQ, but QUAL works here:
                GQ = float(split[5])
            if float(split[5]) > 30 and GQ > 30:
                chrom = split[0]
                locus = int(split[1])
                if ':' not in split[9]:
                    genotype = split[3]
                else:
                    # now for the mindfuck:
                    geno_dat = split[9].split(':')
                    likelihoods = [int(x) for x in geno_dat[1].split(',')]
                    if likelihoods.index(min(likelihoods)) == 0:
                        genotype = split[3]
                    elif likelihoods.index(min(likelihoods)) == 2:
                        if len(split[3]) > len(split[4]):
                            genotype = 'd' + split[3]
                        elif ',' not in split[4]:
                            genotype = split[4]
                        else:
                            genotype = split[4].split(',')[0]
                    elif likelihoods.index(min(likelihoods)) == 5:
                        # if index(min) == 5, must be at least two alt alleles,
                        # therefore must be a comma in field
                        if len(split[3]) > len(split[4].split(',')[1]):
                            genotype = 'd' + split[3]
                        else:
                            genotype = split[4].split(',')[1]
                    elif likelihoods.index(min(likelihoods)) == 9:
                        # as above: must be a comma in field
                        if len(split[3]) > len(split[4].split(',')[2]):
                            genotype = 'd' + split[3]
                        else:
                            genotype = split[4].split(',')[2]


                if chrom in dic:
                    dic[chrom][locus] = genotype
                else:
                    dic[chrom] = {locus: genotype}
    
    n_typed = 0
    for c in d:             # chrom in d
        for l in d[c]:      # locus in d[chrom]
            if c in dic:            # note ref to dic, not d
                if l in dic[c]:     # note ref to dic, not d
                    d[c][l].append(dic[c][l])
                    n_typed += 1
                else:
                    d[c][l].append("N")
            else:
                d[c][l].append("N")

    print f, '\t', n_typed
    return d

def write_all(dic, lis):
    """Write all data from dic to file."""
    names = [x.split('_')[0] for x in lis]
    header = "Locus\t" + '\t'.join(names) + '\n'
    with open(allfil, 'wb') as out:
        out.write(header)
        for chrom in sorted(dic.keys()):
            for loc in sorted(dic[chrom].keys()):
                line = chrom + ':' + str(loc) + '\t' + \
                       '\t'.join(dic[chrom][loc]) + '\n'
                out.write(line)

def write_poly(dic, lis):
    """Write polymorphic data from dic to file."""
    names = [x.split('_')[0] for x in lis]
    header = "Locus\t" + '\t'.join(names) + '\n'
    with open(polyfil, 'wb') as out:
        out.write(header)
        for chrom in sorted(dic.keys()):
            for loc in sorted(dic[chrom].keys()):
                genotypes = set(dic[chrom][loc])
                if 'N' in genotypes:
                    genotypes.remove('N')
                if len(genotypes) > 1:
                    line = chrom + ':' + str(loc) + '\t' + \
                       '\t'.join(dic[chrom][loc]) + '\n'
                    out.write(line)

def write_prop(dic, lis):
    """Write data below propcut from dic to file."""
    cut = round((propcut * 18),0)
    names = [x.split('_')[0] for x in lis]
    header = "Locus\t" + '\t'.join(names) + '\n'
    with open(propfil, 'wb') as out:
        out.write(header)
        for chrom in sorted(dic.keys()):
            for loc in sorted(dic[chrom].keys()):
                genotypes = dic[chrom][loc]
                geno_set = set(genotypes)
                if 'N' in geno_set:
                    geno_set.remove('N')
                if len(geno_set) > 1:
                    if genotypes.count('N') <= (18 - cut):
                        line = chrom + ':' + str(loc) + '\t' + \
                           '\t'.join(dic[chrom][loc]) + '\n'
                        out.write(line)

def write_alltyped(dic, lis):
    """Write locus data if all lines typed from dic to file."""
    names = [x.split('_')[0] for x in lis]
    header = "Locus\t" + '\t'.join(names) + '\n'
    with open(alltyped, 'wb') as out:
        out.write(header)
        for chrom in sorted(dic.keys()):
            for loc in sorted(dic[chrom].keys()):
                genotypes = dic[chrom][loc]
                if genotypes.count('N') == 0:
                    set_g = set(genotypes)
                    if len(set_g) > 1:
                        line = chrom + ':' + str(loc) + '\t' + \
                           '\t'.join(dic[chrom][loc]) + '\n'
                        out.write(line)

if __name__ == '__main__':
    if len(sys.argv) < 6:
        print main.__doc__
        sys.exit()

    allfil = sys.argv[1]
    polyfil = sys.argv[2]
    propcut = float(sys.argv[3])
    propfil = sys.argv[4]
    alltyped = sys.argv[5]
    main()
