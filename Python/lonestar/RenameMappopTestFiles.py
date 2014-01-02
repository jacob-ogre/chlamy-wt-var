# RenameMappopTestFiles.py

import glob
import os

def translation():
	dic = {}
	for line in open("Chlamy_mappop_test_barcodes.tab"):
		if not line.startswith("Barcode"):
			split = line.rstrip().split('\t')
			dic[split[-1]] = split[1]
	return dic


TRANS = translation()
for fil in glob.glob("*.fastq"):
	try:
		split = fil.split('_')
		new = "Chlamy_" + TRANS[split[2]] + ".fastq"
		os.rename(fil, new)
    except:
        print "No TRANS for ", fil