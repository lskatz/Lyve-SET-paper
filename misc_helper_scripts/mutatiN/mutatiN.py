#! /usr/bin/env python

# Author: Amanda Jo Williams-Newkirk, igy7@cdc.gov
# Version: 0.1
# Last updated: 12/17/2015
# Written in Python 2.7.3

# Required packages
import argparse
import os
import random
from itertools import izip as zip, count

# Command line arguments
opt_parser = argparse.ArgumentParser(description = 'Mutate sequences in alignment by inserting Ns at given frequency.')
opt_parser.add_argument('-i', '--infile', dest = 'infile', required = True, help = 'Input file containing aligned sequences in fasta format.')
opt_parser.add_argument('-o', '--outfile', dest = 'outfile', default = os.getcwd() + "/mutated_output", help = 'Output file name.')
opt_parser.add_argument('-f', '--freq', dest = 'freq', default = 0.1, help = 'Float - frequency of Ns in each sequence.')
opt_parser.add_argument('-s', '--seed', dest = 'seed', default = 1, help = 'Integer - set random number generator seed.')
opt_parser.add_argument('-g', '--gaps', dest = 'gaps', default = "F", help = 'Logical (T/F) - preserve existing gap characters in alignment?')
args = opt_parser.parse_args()

# Read in sequences
seqs = dict()
with open(args.infile, 'r') as ifile:
    for line in ifile:
        if ">" in line:
            tmp = line.strip(">")
            tmp = tmp.strip("\n")
            seqs[tmp] = ""
        else:
           seqs[tmp] = seqs[tmp] + line.strip("\n")

# Transform sequences from string to list
for s in seqs:
    seqs[s] = list(seqs[s])

# Mutate sequences
random.seed(int(args.seed))
if args.gaps == "F":
    for s in seqs:
        # Choose locations for Ns
        pos = random.sample(xrange(len(seqs[s])), int(round(float(args.freq) * len(seqs[s]))))
        # Change the positions to Ns
        for p in pos:
            seqs[s][p] = "N"
elif args.gaps == "T":
    for s in seqs:
        # Get non-gap indices for sequence
        notgaps = [i for i, j in zip(count(), seqs[s]) if j != '-']
        # Choose locations for Ns
        pos = random.sample(notgaps, int(round(float(args.freq) * len(notgaps))))
        # Change the positions to Ns
        for p in pos:
            seqs[s][p] = "N"

# Write mutated sequences to output file
with open(args.outfile, 'w') as ofile:
    for s in seqs:
        ofile.write(">" + s + "\n" + "".join(seqs[s]) + "\n")
        