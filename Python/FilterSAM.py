# FilterSAM.py
# Filter SAM output to remove low-quality mappings
# Copyright (C) 2013 Jacob Malcom

# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program; if not, see <http://www.gnu.org/licenses/>.

import glob
from multiprocessing import Pool
import os
import sys

def main():
    """
    Filter SAM output to remove low-quality mappings.

    USAGE:
	python FilterSAM.py <tol> <outdir>

    ARGS:
	tol, the Phred MAPQ above which mappings are retained
	outdir, the directory to which filtered SAM files are written

    RETURNS:
	SAM files from pwd that are filtered to remove poor mappings.

    COMMENTS:
	tol=16 ~ p=0.05
	tol=20 ~ p=0.01
	tol=30 ~ p=0.001
    """
    os.chdir(".")
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    sam_fils = glob.glob("*.sam")

    pool = Pool(processes=12)
    pool.map(filter_sam, sam_fils)

def filter_sam(fil):
    new_fil = outdir + "/filtered_" + fil 
    with open(new_fil, 'wb') as out:
        for line in open(fil):
            if line.startswith('@'):
                out.write(line)
            else:
                if int(line.split('\t')[4]) > tol:
                    out.write(line)

if __name__ == '__main__':
    if len(sys.argv) < 3:
        print main.__doc__
        sys.exit 
    tol = int(sys.argv[1])
    outdir = sys.argv[2]
    main()
