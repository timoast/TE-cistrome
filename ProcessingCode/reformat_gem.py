#! /usr/bin/env python

import sys


with open(sys.argv[1], 'r') as infile:
    header = infile.readline().rsplit()
    print("chr\tstart\tend\t" + "\t".join(header[1:]))
    for line in infile:
        line = line.rsplit()
        coords = line[0].split(":")
        chromosome = "chr" + coords[0]
        position = coords[1]
        end = int(position) + 50
        data = [chromosome, position, str(end)]
        print("\t".join(data) + "\t" + "\t".join(line[1:]))
