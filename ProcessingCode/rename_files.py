#! /usr/bin/env python

import os
from glob import glob

with open("../ProcessedData/accession_names.tsv", "r") as infile:
    for line in infile:
        line = line.rsplit()
        sra = line[0]
        sample = line[1]
        if (os.path.isdir(sra)):
            os.rename(sra, sample)
            os.chdir("./" + sample)
            matching_files = glob(sra + "*")
            for f in matching_files:
                trailing_characters = f.split(sra)[1]
                os.rename(f, sample + trailing_characters)
            os.chdir("..")
        else:
            pass
