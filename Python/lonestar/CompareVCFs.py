# CompareVCFs.py
# Compare two SAMtools .vcf files.
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

import os
from pprint import pprint as pp 
import sys

def main():
    """
    """
    f1_dic = extract_loci()
    # pp(f1_dic)
    compare_fil2(f1_dic)

def extract_loci():
    dic = {}
    for line in open(fil1):
        if not line.startswith("#"):
            split = line.split('\t')
            if split[0] in dic:
                dic[split[0]][split[1]] = {"ref": split[3],
                                           "alt": split[4],
                                           "qual": split[5]}
            else:
                dic[split[0]] = {split[1]: {"ref": split[3],
                                           "alt": split[4],
                                           "qual": split[5]}}
    return dic

def compare_fil2(dic):
    for line in open(fil2):
        if not line.startswith("#"):
            split = line.split('\t')
            if split[1] not in dic[split[0]]:
                print split[0], split[1]
                print split[3], split[4], split[5]

if __name__ == '__main__':
    if len(sys.argv) < 3:
        print main.__doc__
        sys.exit()
    fil1 = sys.argv[1]
    fil2 = sys.argv[2]
    main()