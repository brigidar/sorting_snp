#!/usr/bin/env python


import argparse, os, sys, csv, IPython
import pandas
import pdb
import Bio
from pandas import *
from IPython import get_ipython
from Bio import SeqIO


#output and input file name to give with the script
parser = argparse.ArgumentParser()

parser.add_argument('-o', '--output', help="multifasta with snp as sequence")
parser.add_argument('-s', '--snp_table', help="snp table to transform in fasta")

args = parser.parse_args()
output_file = args.output
input_file = args.snp_table

with open(input_file,'rU') as input:
    with open(output_file,'w') as output:
        sequences = SeqIO.parse(input, "tab")
        count = SeqIO.write(sequences, output, "fasta")

output.close()
input.close()


