#!/usr/bin/python
#
# MakeRADCommandsFile.py

import os
import glob

def main():
    pwd = os.getcwd()
    gmap_cmd = pwd + "/RADgmapcommands"     # gmapper cmds
    dat_dir = "/scratch/01703/jmalcom/Chlamy_wt_data/"
    sam_dir = "/scratch/01703/jmalcom/Chlamy_RAD_SAM/"

    with open(gmap_cmd, 'wb') as o:
        for fil in glob.glob(dat_dir + "filt_cat*.fastq"):
            fil_only = fil.split('/')[-1]
            sample = fil_only.split('_')[-1].split('.')[0]
            gmap = "Ch_" + sample + "_fq.sam"
            
            # make the line for calling gmapper-cs and write
            line = "$SH23DIR/gmapper-cs -N 24 -o 10 -Q --qv-offset 33 $DATADIR/" \
+ fil_only + " -L $REFPROJ/Chlamy_genome_projection_cs > $GMAPDIR/" + gmap + ";\n"

            o.write(line)

        

if __name__ == '__main__':
    main()
