# CompileAllGenotypes.py
# Compile the .vcf files into a single table.
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
    dic1 = {}

    for line in open("Chlamy_wt_poly_loci.tab"):
        split = line.split('\t')
        dic1[split[0]] = split[1:]

    for line in open("Chlamy_wt_poly_loci_v2.tab"):
        split = line.split('\t')
        if split[0] not in dic1:
            print line

if __name__ == '__main__':
    
    main()
