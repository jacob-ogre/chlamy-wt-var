#!/usr/bin/python
#
# MakeRADCommandsFile.py

import os
import glob

def main():
    pwd = os.getcwd()
    sam_cmd = pwd + "/SAMcommands"     # gmapper cmds
    sam_dir = "/scratch/01703/jmalcom/Chlamy_RAD_SAM/"

    list_lines = get_lines(sam_dir)
    os.chdir(pwd)
    with open(sam_cmd, 'wb') as out:
        write_view(list_lines, out)
        write_sort(list_lines, out)
        write_index(list_lines, out)
        write_mpileup(list_lines, out)
        write_bcftool(list_lines, out)

def get_lines(sdir):
    list_lines = []
    os.chdir(sdir)
    for fil in glob.glob("*.sam"):
        list_lines.append(fil.split('_')[1])
    return list_lines

def write_view(lis, out):
    for i in lis:
        line = "$SAMTL/samtools view -b -q 20 -S $CRDAT/Ch_" + i + "_fq.sam > $CRDAT/" + i + "_filtered.bam;\n"
        out.write(line)
        out.write("wait;\n")

def write_sort(lis, out):
    for i in lis:
        line = "$SAMTL/samtools sort $CRDAT/" + i + "_filtered.bam $CRDAT/" + i +"_filtered_sorted;\n"
        out.write(line)
        out.write("wait;\n")

def write_index(lis, out):
    for i in lis:
        line = "$SAMTL/samtools index $CRDAT/" + i + "_filtered_sorted.bam;\n"
        out.write(line)
        out.write("wait;\n")

def write_mpileup(lis, out):
    for i in lis:
        line = "$SAMTL/samtools mpileup -uEf $CRREF/Cre_comp_ref.fa $CRDAT/" + i + "_filtered_sorted.bam > $CRDAT/" + i + ".bcf;\n"
        out.write(line)
        out.write("wait;\n")

def write_bcftool(lis, out):
    for i in lis:
        line = "$SAMTL/bcftools/bcftools view -cg $CRDAT/" + i + ".bcf > $CRDAT/" + i + ".vcf;\n"
        out.write(line)
        out.write("wait;\n")

if __name__ == '__main__':
    main()
