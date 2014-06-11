# split_BLUPs_file.py
# Convert BLUPs file to a line x environment matrix.
# Copyright (C) 2014 Jacob Malcom, jmalcom@uconn.edu
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.

import sys

def main():
    """Convert BLUPs file to a line x environment matrix.

    USAGE:
        python split_BLUPs_file.py <infile> <outfil>
    ARGS:
        infile, path to a tab'd file of BLUPs from R
        outfil, path to the tab'd file to be written
    """
    environs = ["line"]
    line_dat = {}
    for line in open(infile):
        if not line.startswith("(Intercept)"):
            data = line.rstrip().split("\t")
            first = data[0].split(":")
            if first[0] in line_dat:
                line_dat[first[0]].append(data[1])
            else:
                line_dat[first[0]] = [data[1]]
            if first[1] not in environs:
                environs.append(first[1])

    with open(outfil, 'w') as out:
        header = "\t".join(environs)
        out.write(header + "\n")
        for k in sorted(line_dat):
            line_dat[k].insert(0, k)
            to_write = "\t".join(line_dat[k]) + "\n"
            out.write(to_write)

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print main.__doc__
        sys.exit()
    infile = sys.argv[1]
    outfil = sys.argv[2]
    main()