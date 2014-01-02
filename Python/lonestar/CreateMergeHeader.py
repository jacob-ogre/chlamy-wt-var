# CreateMergeHeader.py
# Creates a .sam header for merging two or more files with samtools merge.
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
# import os
import sys

def main():
    """
    Creates a .sam header for merging two or more files with samtools merge.

    USAGE:
        python CreateMergeHeader.py <filebase> <outdir>

    ARGS:
        filebase, the base name of the sam/bam files to be merged
        outdir, the directory into which the updated header should be written

    RETURNS:
        outdir/<filebase>_header.sam 

    COMMENTS:
            *********************** NOTE ***********************
            This program assumes that one short-read mapper was 
            used to create the sam files; if that is not the 
            case, then the user needs to edit lines 92-96 to 
            accomodate multiple @PG flags.
            ****************************************************

        GATK expects @RG data for every read group present in a file, so using 
        SAMtools merge without the -h flag--which indicates the header file to 
        use to put at the top of the merged file--will fail. This script creates 
        a merged header file in the same directory as the .sam files.
    """
    pattern = filebase + "*.sam"
    header_dat = {}
    for fil in glob.glob(pattern):
        header_dat = extract_header_info(fil, header_dat)
    write_header(header_dat)

def extract_header_info(fil, dic):
    """
    Store each line from the header in dic as a list element.

    COMMENTS:
        Most lines (e.g., @SQ) are idential between .sam files (or they should 
        be if you've done everything correctly!), but the @RG and @PG lines
        differ. BREAK to avoid processing the entire .sam.
    """
    line_numb = 0
    for line in open(fil):
        if line.startswith("@"):
            if line_numb in dic:
                dic[line_numb].append(line)
            else:
                dic[line_numb] = [line]
            line_numb += 1
        else:
            break
    return dic

def write_header(dic):
    """
    Write a single header that includes @RG for each sample.

    COMMENTS:
        The key is to write one copy of the lines that are the same between 
        files (e.g. @SQ lines) but to then write each @RG or @PG line so that
        GATK has all of the information it needs. Test the need for only one
        copy by calculating the set of the list of lines in dic.
    """
    outfile = outdir + '/' + filebase + "_header.sam"
    with open(outfile, 'wb') as out:
        for line in sorted(dic.keys()):
            if len(set(dic[line])) == 1:
                out.write(dic[line][0])
            elif dic[line][0].startswith("@RG"):
                for i in dic[line]:
                    out.write(i)
            else:
                out.write(dic[line][0])

if __name__ == '__main__':
    if len(sys.argv) < 3:
        print main.__doc__
        sys.exit()

    filebase = sys.argv[1]
    outdir = sys.argv[2]
    main()