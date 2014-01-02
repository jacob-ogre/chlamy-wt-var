# MakeSAM_GATKCommands.py
# Generates commands for running both SAM and GATK on a set of .bam files.
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
import sys

def main():
    """
    Create a commands file for SAM and GATK genotyping using TACC Launcher.

    USAGE:
        python MakeSAM_GATKCommands.py <srcdir> <outfil>

    ARGS:
        srcdir, the path to the .sam files to be processed
        outfil, the name of the launcher commands file

    RETURNS:
        outfil, which is used by the launcher .sge file

    COMMENTS:
        Globs for .sam files using the patterns *.sam.  Throughout, this 
        code-generating script assumes that the sample identifier is enclosed 
        between two '_' characters, next-to-last in the file path (see line 45).
        Note that the path prefixes are hard-coded here: make adjustments in the
        .sge file that uses TACC Launcher.
    """
    with open(outfil, "wb") as out:
        counter = 0
        for fil in glob.glob(srcdir + "*.sam"):
            base = fil.split('.')[0].split('/')[-1]
            # line = base.split('_')[0]
            # lane = base.split('_')[-1]

            # # First, the picard RG additions:
            # cd1 = "cd $PICARD && pwd && "
            # out.write(cd1)

            # piRG = "java -Xmx1G -jar AddOrReplaceReadGroups.jar INPUT=" + fil + \
            #        " OUTPUT=$UPDAT/" + base + "_updated.sam " + \
            #        "RGID=\"" + str(counter) + "\" " + \
            #        "RGLB=\"Aug12\" " + \
            #        "RGPL=\"solid\" " + \
            #        "RGPU=\"" + lane + "\" " + \
            #        "RGSM=\"" + line + "\" && "
            # out.write(piRG)

            # # Second, the samtools additions:
            cd2 = "cd $SAMTL && "
            out.write(cd2)

            # stVW = "samtools view -b -q 20 -S $UPDAT/" + base + \
            #        "_updated.sam  " + \
            #        "> $BAMS/" + base + ".bam && "
            # out.write(stVW)

            # stSO = "samtools sort $BAMS/" + base + ".bam $SORT/" + base + \
            #        "_sorted \n"
            # out.write(stSO)

            base = base.split('_')[0]

            stIN = "samtools index $MERGE/" + base + ".bam && "
            out.write(stIN)

            stMP = "samtools mpileup -uEf $CRREF/Cre_comp_ref.fa " + \
                   "$MERGE/" + base + ".bam > $BCFS/" + base + ".bcf && "
            out.write(stMP)

            stBC = "bcftools/bcftools view -cg $BCFS/" + base + ".bcf > " + \
                   "$VCFS/" + base + "_SAM.vcf && "
            out.write(stBC)

            # Third, GATK analysis:
            cd3 = "cd $WORK/GATK && pwd && "
            out.write(cd3)

            # realign around indels:
            gkRC = "java -Xmx1G -jar GenomeAnalysisTK.jar -T RealignerTargetCreator  " + \
                   "-I $MERGE/" + base + ".bam  " + \
                   "-R $CRREF/Cre_comp_ref.fa  " + \
                   "-o $REAL/" + base + "_realigner.intervals && "
            out.write(gkRC)

            gkIR = "java -Xmx1G -jar GenomeAnalysisTK.jar -T IndelRealigner " + \
                   "-I $MERGE/" + base + ".bam " + \
                   "-R $CRREF/Cre_comp_ref.fa " + \
                   "-targetIntervals $REAL/" + base + "_realigner.intervals " + \
                   "-o $RBAM/" + base + "_realigned.bam && "
            out.write(gkIR)

            # first pass of the Unified Genotyper to get high-quality variants
            gkUG = "java -Xmx1G -jar GenomeAnalysisTK.jar " + \
                   "-T UnifiedGenotyper -nt 12 " + \
                   "-R $CRREF/Cre_comp_ref.fa  " + \
                   "-I $RBAM/" + base + "_realigned.bam  " + \
                   "-ploidy $PLOIDY " + \
                   "-glm BOTH " + \
                   "-out_mode EMIT_ALL_CONFIDENT_SITES " + \
                   "-o $VCFS/" + base + "_GATK.vcf && "
            out.write(gkUG)

            # Filter out all but the highest-quality SNPs given GQ > 90 and 
            # MAPQ > 90th percentile for the sample:
            pyHQ = "python $VCFS/VCFHighQualityVariants.py " + \
                   "$VCFS/" + base + "_GATK.vcf " + \
                   "$SNPS && "
            out.write(pyHQ)

            # Now do the base quality score recalibration for the sample:
            gkBR = "java -Xmx1G -jar GenomeAnalysisTK.jar -T BaseRecalibrator " + \
                   "-I $RBAM/" + base + "_realigned.bam " + \
                   "-R $CRREF/Cre_comp_ref.fa " + \
                   "-knownSites $SNPS/" + base + "_HQ_SNPs.vcf " + \
                   "-o $RECAL/" + base + "_recal_data.grp && "
            out.write(gkBR)

            gkPR = "java -Xmx1G -jar GenomeAnalysisTK.jar -T PrintReads " + \
                   "-I $RBAM/" + base + "_realigned.bam " + \
                   "-R $CRREF/Cre_comp_ref.fa " + \
                   "-BQSR $RECAL/" + base + "_recal_data.grp " + \
                   "-o $REREBAM/" + base + "_realigned_recal.bam && "
            out.write(gkPR)

            # And the final GATK UnifiedGenotyper run with corrected quals:
            gkG2 = "java -Xmx1G -jar GenomeAnalysisTK.jar " + \
                   "-T UnifiedGenotyper -nt 12 " + \
                   "-R $CRREF/Cre_comp_ref.fa " + \
                   "-I $REREBAM/" + base + "_realigned_recal.bam " \
                   "-ploidy $PLOIDY " + \
                   "-glm BOTH " + \
                   "-out_mode EMIT_ALL_CONFIDENT_SITES " + \
                   "-o $VCFS/" + base + "_GATK_recal.vcf \n"
            out.write(gkG2)

            counter += 1



if __name__ == '__main__':
    if len(sys.argv) < 3:
        print main.__doc__
        sys.exit()

    srcdir = sys.argv[1]
    outfil = sys.argv[2]
    main()