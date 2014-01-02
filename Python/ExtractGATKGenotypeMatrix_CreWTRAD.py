# ExtractGATKGenotypeMatrix.py
# Compile GATK .vcf files into a single table.
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
                            @@@@@ NOTE! @@@@@@
            This script is only appropriate for non-merged sam/bam files! If you
            have merged sam/bam during the genotyping workflow, please edit the 
            split[9] != '.' to accomodate as many extra columns as necessary. 
            For example, if you merged two bams with samtools merge, then called
            genotypes with GATK, you will need have calls for split[9] and 
            split[10].
                            @@@@@@@@@@@@@@@@@

        Uses glob on the pattern *filtered.vcf to pull in .vcf files; note that
        the filtered.vcf files have already had low-quality genotype calls 
        removed.  Although each of the write_* functions could be cleaned up to 
        shorten the script, it's still pretty fast (~5 min currently).  
    """
    os.chdir(os.getcwd())

    all_loci = {}
    fil_list = glob.glob(indir + "*GATK_recal.vcf")
    for fil in fil_list:
        all_loci = get_all_loci(fil, all_loci)

    for fil in fil_list:
        all_loci = fill_dict(fil, all_loci)

    write_all(all_loci, fil_list)

def get_all_loci(f, d):
    """Read loci from file f and add to dict d."""
    for line in open(f):
        if not line.startswith('#'):
            split = line.split('\t')
            if split[9] != '.':
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
            chrom = split[0]
            locus = split[1]
            if split[4] == '.':
                if len(split[3]) > 1:
                    locus = locus + 'd'

                if ':' in split[-1]:
                    depth = int(split[-1].split(':')[1])
                else:
                    depth = 1

                if depth > 2 and len(split[3]) == 1:
                    genotype = split[3]
                elif depth > 2 and len(split[3]) > 1:
                    genotype = 'd' + split[3]
                else:
                    genotype = "N"

            else:
                # we prepend a 'd' to indicate a deletion:
                if len(split[3]) > 1:
                    genotype = 'd' + split[3]
                    locus = locus + 'd'

                # rather than throw an exception, just print if >1 ALT allele.
                elif ',' in split[4]:
                    print f, '\n', line

                else:
                    if ':' in split[-1]:
                        cur_dat = split[-1].split(':')
                        if len(cur_dat) == 2:
                            depth = int(cur_dat[1])
                        else:
                            depth = int(cur_dat[2])
                    else:
                        depth = 1

                    if depth > 2 and len(split[4]) > 1 and ',' not in split[4]:
                        genotype = 'i' + split[4]
                        locus = locus + 'i'
                    elif depth > 2 and len(split[4]) == 1:
                        genotype = split[4]
                    else:
                        genotype = "N"

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
    header = "Scaffold:Locus\t" + '\t'.join(names) + '\n'
    with open(allfil, 'wb') as out:
        out.write(header)
        for chrom in sorted(dic.keys()):
            for loc in sorted(dic[chrom].keys()):
                if dic[chrom][loc].count("N") != len(dic[chrom][loc]):
                    line = chrom + ':' + str(loc) + '\t' + \
                           '\t'.join(dic[chrom][loc]) + '\n'
                    out.write(line)


if __name__ == '__main__':
    if len(sys.argv) < 3:
        print main.__doc__
        sys.exit()

    indir = sys.argv[1]
    allfil = sys.argv[2]
    main()
