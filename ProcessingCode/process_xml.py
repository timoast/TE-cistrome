#! /usr/bin/env python

"""
Extract the sample names from study xml file, pair with SRA accession names
to allow meaningful names to be assigned to the downloaded sequencing data
"""

import xml.etree.cElementTree as ET


tree = ET.parse('../RawData/SraExperimentPackage.xml')
root = tree.getroot()

title = root.getiterator("TITLE")
runs = root.getiterator("RUN")

sample_names = [i.text for i in title]
accession_numbers = [i.attrib["accession"] for i in runs]

for i in sample_names:
    if i.startswith("GSM"):
        sample_names.remove(i)

with open("../ProcessedData/accession_names.tsv", "w+") as outfile:
    for i in xrange(len(sample_names)):
        outfile.write(accession_numbers[i] + "\t" + sample_names[i] + "\n")
