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
#

import glob
import os
import sys

def main():
    """
    Compile genotypes from multiple .vcf files into line x genotype .tab files.

    USAGE:
        python CompileAllGenotypes.py <inpath> <allfil>
    """
    os.chdir(os.getcwd())

    all_loci = {}
    fil_list = glob.glob(inpath + "*GATK_recal.vcf")
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
                locus = split[1]
                if len(split[3]) > 1:
                    locus = locus + "d"
                if len(split[4]) > 1 and "," not in split[4]:
                    locus = locus + "i"
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

                if depth > minD and len(split[3]) == 1:
                    genotype = split[3]
                elif depth > minD and len(split[3]) > 1:
                    genotype = 'd' + split[3]
                else:
                    genotype = "N"

            else:
                if len(split[3]) > 1:
                    genotype = 'd' + split[3]
                    locus = locus + 'd'
                elif ',' in split[4]:
                    print f, '\n', line
                    continue
                else:
                    qual = 20
                    if ':' in split[-1]:
                        cur_dat = split[-1].split(':')
                        if len(cur_dat) == 2 or len(cur_dat) == 3:
                            depth = int(cur_dat[1])
                        elif len(cur_dat) > 3: 
                            al_counts = [int(x) for x in cur_dat[1].split(',')]
                            try:
                                if min(al_counts) / sum(al_counts) > 0.2:
                                    continue
                            except:
                                continue
                            qual = int(cur_dat[3])
                    else:
                        depth = 1

                    if depth > minD and qual > minQ and len(split[4]) > 1:
                        genotype = "i" + split[4]
                        locus = locus + "i"
                    elif depth > minD and qual > minQ and len(split[4]) == 1:
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
    """Write locus data if all lines typed from dic to file."""
    fil_list = [x.split("/")[-1] for x in lis]
    names = [x.split('_')[0] for x in fil_list]
    header = "Scaffold:Locus\t" + '\t'.join(names) + '\n'
    with open(allfil, 'wb') as out:
        out.write(header)
        for chrom in sorted(dic.keys()):
            for loc in sorted(dic[chrom].keys()):
                genotypes = dic[chrom][loc]
                if genotypes.count("N") < 16:
                    line = chrom + ':' + loc + '\t' + \
                       '\t'.join(genotypes) + '\n'
                    out.write(line)

if __name__ == '__main__':
    if len(sys.argv) != 5:
        print main.__doc__
        sys.exit()

    minD = int(sys.argv[1])
    minQ = int(sys.argv[2])
    inpath = sys.argv[3]
    allfil = sys.argv[4]
    main()
