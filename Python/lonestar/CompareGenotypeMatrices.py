# CompareGenotypeMatrices.py
# Determine the degree to which two genotype matrices are similar.
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

# import os
import sys

def main():
    """
    """
    fil1_dict, lines = extract_data(fil1)
    fil2_dict, lines = extract_data(fil2)
    compare_dicts(fil1_dict, fil2_dict, lines)

def extract_data(fil):
    dic = {}
    for line in open(fil):
        if not line.startswith("Loc"):
            split = line.rstrip().split('\t')
            dic[split[0]] = split[1:]
        else:
            lis = line.rstrip().split('\t')[1:]
    return dic, lis

def compare_dicts(d1, d2, lis):
    n_same = 0
    n_diff = 0
    both_n = 0
    one_typed = 0
    n_missing_2 = 0
    n_missing_1 = 0
    for k in d1:
        if k in d2:
            if d1[k] == d2[k]:
                n_same += 18
            else:
                for i in range(len(d1[k])):
                    if d1[k][i] == "N" and d2[k][i] == "N":
                        both_n += 1
                    elif d1[k][i] == d2[k][i]:
                        n_same += 1
                    elif d1[k][i] == "N" or d2[k][i] == "N":
                        one_typed += 1
                    else:
                        n_diff += 1
                        print k, lis[i], d1[k][i], d2[k][i]
        else:
            n_missing_2 += 1

    for k in d2:
        if k not in d1:
            n_missing_1 += 1

    print fil1, fil2
    print "# same:", n_same
    print "# diff:", n_diff
    print "# both N:", both_n
    print "# only one typed:", one_typed
    print "# loci in", fil1, "not in", fil2, n_missing_2
    print "# loci in", fil2, "not in", fil1, n_missing_1


if __name__ == '__main__':
    if len(sys.argv) < 3:
        print main.__doc__
        sys.exit()

    fil1 = sys.argv[1]
    fil2 = sys.argv[2]
    main()
