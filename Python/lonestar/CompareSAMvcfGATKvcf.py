# CompareSAMvcfGATKvcf.py
# Compare the 2334 SAM-generated vcf to the GATK-generated vcf
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

number_shared = 0
number_diff = 0
inSAMnotGATK = 0
inGATKnotSAM = 0

dic = {}
for line in open("vcfs/filtered_vcfs/2343_SAM_filtered.vcf"):
    if not line.startswith('#'):
        split = line.split('\t')
        chrom = split[0]
        locus = int(split[1])

        # if REF, else ALT:
        if ':' not in split[9]:
            genotype = split[3]
        else:
            geno_dat = split[9].split(':')
            likelihoods = [int(x) for x in geno_dat[1].split(',')]
            if likelihoods.index(min(likelihoods)) == 0:
                # this shouldn't be necessary, but just in case
                genotype = split[3]
            elif likelihoods.index(min(likelihoods)) == 2:
                if ',' not in split[4]:
                    genotype = split[4]
                else:
                    genotype = split[4].split(',')[0]
            elif likelihoods.index(min(likelihoods)) == 5:
                # if index(min) == 5, must be at least two alt alleles,
                # therefore must be a comma in field
                genotype = split[4].split(',')[1]
            elif likelihoods.index(min(likelihoods)) == 9:
                # as above: must be a comma in field
                genotype = split[4].split(',')[2]
            else:
                # a catch-all in case I missed something:
                genotype = "N"
                print line

        if chrom in dic:
            dic[chrom][locus] = genotype
        else:
            dic[chrom] = {locus: genotype}

GATKdic = {}
for line in open("2343_GATK_recal.vcf"):
    if not line.startswith('#'):
        split = line.rstrip().split('\t')
        chrom = split[0]
        locus = int(split[1])
        if split[9] != '.':
            gen_dat = split[9].split(':')[0]

            if gen_dat == '0':
                genotype = split[3]
            else:
                genotype = split[4]

            if chrom in GATKdic:
                GATKdic[chrom][locus] = genotype
            else:
                GATKdic[chrom] = {locus: genotype}

with open("InSAMnotGATK.out", 'wb') as out:
    for chrom in dic:
        for locus in dic[chrom]:
            if chrom in GATKdic:
                if locus in GATKdic[chrom]:
                    if dic[chrom][locus] == GATKdic[chrom][locus]:
                        number_shared += 1
                    else:
                        number_diff += 1
                        print chrom, locus, dic[chrom][locus], GATKdic[chrom][locus]
                else:
                    inSAMnotGATK += 1
                    line = chrom + '\t' + str(locus) + '\n'
                    out.write(line)
            else:
                inSAMnotGATK += 1
                line = chrom + '\t' + str(locus) + '\n'
                out.write(line)
    print "#################"
    print "#################"
    print "#################"
    print "#################"
    print "#################"

with open("InGATKnotSAM.out", 'wb') as out:
    for chrom in GATKdic:
        for locus in GATKdic[chrom]:
            if chrom in dic:
                if locus not in dic[chrom]:
                    number_diff += 1
                    inGATKnotSAM += 1
                    line = chrom + '\t' + str(locus) + '\n'
                    out.write(line)
            else:
                number_diff += 1
                inGATKnotSAM += 1
                line = chrom + '\t' + str(locus) + '\n'
                out.write(line)


print "N_shared", number_shared 
print "N_diff", number_diff
print "N_inSAMnotGATK", inSAMnotGATK
print "N_inGATKnotSAM", inGATKnotSAM