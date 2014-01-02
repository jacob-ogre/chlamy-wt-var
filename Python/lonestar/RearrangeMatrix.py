# RearrangeMatrix.py
# One-off script to rearrange in-house genotype matrix
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

import sys

def main():
    """USAGE: python RearrangeMatrix.py <infile> <outfile>"""
    with open(outfile, "wb") as out:
        for line in open(infile):
            if not line.startswith("Loci"):
                split = line.rstrip().split('\t')
                move = split.pop(16)
                split.append(move)
                new_line = '\t'.join(split) + '\n'
                out.write(new_line)
            else:
                split = line.rstrip().split('\t')
                split.pop(16)
                split.append("CHL90")
                new_line = '\t'.join(split) + '\n'
                out.write(new_line)


if __name__ == '__main__':
    if len(sys.argv) < 3:
        print main.__doc__
        sys.exit

    infile = sys.argv[1]
    outfile = sys.argv[2]
    main()