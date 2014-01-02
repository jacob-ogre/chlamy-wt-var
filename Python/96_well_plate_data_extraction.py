# 96_well_plate_extraction_v3.py
# New version of 96-well plate extration.
# Copyright (C) 2013 Jacob Malcom, jmalcom@uconn.edu

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
from pprint import pprint as pp
import sys

ROWS = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']
COLS = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12']
STDS = ["REF0", "REF1", "REF10", "REF20", "REF5", "REF50", "N++", "TAP"]
LOCATIONS = {"B2": 1, "B3": 2, "B4": 3, "B5": 4, "B6": 5,
             "C2": 6, "C3": 7, "C4": 8, "C5": 9, "C6": 10,
             "D2": 11, "D3": 12, "D4": 13, "D5": 14, "D6": 15,
             "G2": 1, "G3": 2, "G4": 3, "G5": 4, "G6": 5,
             "F2": 6, "F3": 7, "F4": 8, "F5": 9, "F6": 10,
             "E2": 11, "E3": 12, "E4": 13, "E5": 14, "E6": 15,
             "B11": 1, "B10": 2, "B9": 3, "B8": 4, "B7": 5,
             "C11": 6, "C10": 7, "C9": 8, "C8": 9, "C7": 10,
             "D11": 11, "D10": 12, "D9": 13, "D8": 14, "D7": 15,
             "G11": 1, "G10": 2, "G9": 3, "G8": 4, "G7": 5,
             "F11": 6, "F10": 7, "F9": 8, "F8": 9, "F7": 10,
             "E11": 11, "E10": 12, "E9": 13, "E8": 14, "E7": 15}


def main():
    """Extract data from SpectraMax plate reader output, write tab-delim. out.

    USAGE
        python 96_well_plate_extraction_v3.py <key_dir> <in_dir> <out_fil>

    ARGS
        key_dir, path to the directory with the plate-line mappings
        in_dir, path to the directory containing raw data
        out_fil, path of the directory for output files
    RETURNS
        A tab-delimited file of plate, line, date, and density data.
    """
    mappings = load_plate_mappings()
    # pp(mappings)
    raw_data = load_plate_data()
    write_data(mappings, raw_data)

def load_plate_mappings():
    """Return dict of form plate:row:col:line/treatment/coord."""
    res = {}
    for fil in glob.glob(key_dir + "*.tab"):
        cur_plate = ""
        for line in open(fil):
            data = line.rstrip().split("\t")
            if not data[:2] == ["", ""] and len(data) > 3:
                if data[1].startswith("PLATE"):
                    cur_plate = data[1][6:]
                elif data[1] in ROWS[1:7]:
                    for i in range(3, len(data)):
                        if cur_plate not in res:
                            res[cur_plate] = {data[1]: {COLS[i-2]:
                                             {"line": strip(data[i])}}}
                        elif data[1] not in res[cur_plate]:
                            res[cur_plate][data[1]] = {COLS[i-2]:
                                                      {"line": strip(data[i])}}
                        elif COLS[i-2] not in res[cur_plate][data[1]]:
                            res[cur_plate][data[1]][COLS[i-2]] = {"line":
                                                                  strip(data[i])}
                        else:
                            res[cur_plate][data[1]][COLS[i-2]]["treat"] = data[i]
                            coord = data[1] + COLS[i-2]
                            res[cur_plate][data[1]][COLS[i-2]]["coord"] = \
                                                       LOCATIONS[coord]

    return res

def strip(x):
    return x.split(".")[0]

def load_plate_data():
    """Return a dict of form plate:row:col:date:data."""
    res = {}
    for fil in glob.glob(in_dir + "*Oct*" + "/*.txt"):
        res = extract_file(fil, res)
    return res

def extract_file(f, d):
    """Return an updated dict d of form plate:row:col:date:data."""
    plate, date = extract_file_data(f)
    line_ct = 0
    if plate not in d:
        d[plate] = {}
    for line in open(f, 'rU'):
        line_ct += 1
        if line_ct > 4 and line_ct < 11:
            cur_row = ROWS[line_ct - 4]
            if cur_row not in d[plate]:
                d[plate][cur_row] = {}
            data = line.rstrip().split('\t')[3:-1]
            for i in range(len(data)):
                if COLS[i+1] not in d[plate][cur_row]:
                    d[plate][cur_row][COLS[i+1]] = {date: data[i]}
                else:
                    d[plate][cur_row][COLS[i+1]][date] = data[i]
    return d

def extract_file_data(f):
    """Return a dict of form "plate": plate, "date": date."""
    split1 = f.split('/')[-1].split('.')[0].split('_')
    d = split1[1].split("-")
    date = d[0] + ' ' + d[1] + ' ' + d[2]
    return split1[0], date

def write_data(m, dat):
    """Write tab-delimited file of extracted plate data."""
    with open(all_fil, 'wb') as out:
        header = "Date\tPlate\tRow\tColumn\tCoord\tEnv\tLine\tFluor\n"
        out.write(header)
        for p in sorted(dat.keys()):
            for r in sorted(dat[p].keys()):
                for c in sorted(dat[p][r].keys()):
                    for d in sorted(dat[p][r][c].keys()):
                        coord = str(m[p][r][c]["coord"])
                        genot = m[p][r][c]["line"]
                        treat = m[p][r][c]["treat"]
                        fitn = dat[p][r][c][d]
                        lin_dat = [d, p, r, c, coord, treat, genot, fitn]
                        out.write(('\t'.join(lin_dat) + '\n'))

    with open(max_fil, 'wb') as out:
        header = "Date\tPlate\tRow\tColumn\tCoord\tEnv\tLine\tFluor\n"
        out.write(header)
        for p in sorted(dat.keys()):
            for r in sorted(dat[p].keys()):
                for c in sorted(dat[p][r].keys()):
                    if m[p][r][c] not in STDS:
                        cur_max = 0
                        max_dat = ""
                        for d in dat[p][r][c]:
                            if float(dat[p][r][c][d]) > cur_max:
                                cur_max = float(dat[p][r][c][d])
                                max_dat = d
                        coord = str(m[p][r][c]["coord"])
                        genot = m[p][r][c]["line"]
                        treat = m[p][r][c]["treat"]
                        lin_dat = [max_dat, p, r, c, coord, treat, genot, str(cur_max)]
                        out.write(('\t'.join(lin_dat) + '\n'))

    with open(tot_fil, 'wb') as out:
        header = "Date\tPlate\tRow\tColumn\tCoord\tEnv\tLine\tTotal\n"
        out.write(header)
        for p in sorted(dat.keys()):
            for r in sorted(dat[p].keys()):
                for c in sorted(dat[p][r].keys()):
                    if m[p][r][c] not in STDS:
                        cur_tot = 0
                        for d in dat[p][r][c]:
                            cur_tot += float(dat[p][r][c][d])
                        coord = str(m[p][r][c]["coord"])
                        genot = m[p][r][c]["line"]
                        treat = m[p][r][c]["treat"]
                        lin_dat = [d, p, r, c, coord, treat, genot, str(cur_tot)]
                        out.write(('\t'.join(lin_dat) + '\n'))


if __name__ == '__main__':
    base = "/Users/jacobmalcom/"
    key_dir = base + "chlamy-mapping-population/Python/wt_fit_plate_layouts_files/"
    in_dir = base + "Dropbox/Chlamy_project/Chlamy_wt_fitness_assays/"
    all_fil = in_dir + "extracted_data/wt_fitness_data_FINAL.tab"
    max_fil = in_dir + "extracted_data/wt_fitness_maxF_FINAL.tab"
    tot_fil = in_dir + "extracted_data/wt_fitness_totF_FINAL.tab"
    main()
