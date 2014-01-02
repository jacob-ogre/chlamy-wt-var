# ConcatenateFilteredReadFiles.py
# Concatenates the many csfasta files into one file per line, for mapping.
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
import shutil

def main():
    """
    Concatenates the many csfasta files into one file per line, for mapping.

    USAGE:
        python ConcatenateFilteredReadFiles.py

    ARGS:
        none explicit; glob used extensively

    RETURNS:
        one file each .csfasta and .qual with prefix filt_cat_Ch*
    """
    os.chdir("/scratch/01703/jmalcom/Chlamy_wt_data/")
    bases = get_base_names()
    for base in bases:
        order = []
        dest_fil = "filt_cat_Ch_" + base + ".csfasta"
        with open(dest_fil, 'wb') as out:
            search = "filt_Ch_" + base + "*.csfasta"
            for fil in glob.glob(search):
                order.append(fil)
                print fil
                shutil.copyfileobj(open(fil, 'rb'), out)
        print order

        dest2 = "filt_cat_Ch_" + base + "_QV.qual"
        with open(dest2, 'wb') as out2:
            for fil in order:
                qual_fil = fil.split('.')[0] + "_QV.qual"
                print qual_fil
                shutil.copyfileobj(open(qual_fil, 'rb'), out2)


def get_base_names():
    list_lines = []
    for fil in glob.glob("*.csfasta"):
        split = fil.split('_')[-1].split('.')[0]
        if len(split) == 3:
            list_lines.append(split[:2])
        else:
            list_lines.append(split[:4])
    return list(set(list_lines))

if __name__ == '__main__':
    main()
