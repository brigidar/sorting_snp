#!/usr/bin/env python



import argparse, os, sys, csv, IPython
import pandas
import pdb
import glob
from pandas import *
from IPython import get_ipython
from glob import glob

#output file name to give with the script
parser = argparse.ArgumentParser()

parser.add_argument('-o', '--output', help="combined snp table and add csv extension to name")

args = parser.parse_args()
output_file = args.output


#calls up all txt files in the current directory

files = glob('*.txt')

# concatenates them all together and reads molecule and refpos as indexes
data = [read_csv(f, sep ='\t', header=0, index_col=[0,1], dtype=unicode) for f in files]
data = [d.T for d in data]
d = concat(data, keys=files)
# transpose the table to have the correct headers
df = d.T


# find columns that have qbase in it and add them together with the position and the refbase and molecule name
count_qbase=df.columns.values
qindexes=[]
for i, v in enumerate(count_qbase):
    if 'qbase:' in v[1]:
     qindexes.append(i)
df2=(df.iloc[:,qindexes]).T
df3=(df.iloc[:,0:2]).T
final=(concat([df3,df2]).T)



# save new csv table

with open(output_file,'w') as output:
 final.to_csv(output)
