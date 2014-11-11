#for replacement of a given value with NaN
#http://stackoverflow.com/questions/18172851/deleting-dataframe-row-in-pandas-based-on-column-value

# to remove any symbol and replace it with nan
#http://stackoverflow.com/questions/875968/how-to-remove-symbols-from-a-string-with-python


#!/usr/bin/env python

#before reading in text file from Galaxy add a row with -- to have a second header

import argparse, os, sys, csv, IPython
import pandas
import pdb
from pandas import *
from IPython import get_ipython
import matplotlib.pyplot as plt
from pandas.util.testing import assert_frame_equal


#output and input file name to give with the script
parser = argparse.ArgumentParser()

parser.add_argument('-o', '--output', help="sorted snp table with identical and no hit removed for snp fasta")
parser.add_argument('-s', '--snp_table', help="snp table to sort")
parser.add_argument('-t', '--total', help="inverted table to look at")


args = parser.parse_args()
output_file = args.output
input_file = args.snp_table
output2_file = args.total


#read in file as dataframe
df =read_csv(input_file, sep='\t', index_col=[0,1], header=0, dtype=unicode)
df=df.drop('syn?', axis=1)

count_qbase=list(df.columns.values)
qindexes=[]
for i, v in enumerate(count_qbase):
    if 'qbase:' in v:
        qindexes.append(i)
df2=(df.iloc[:,qindexes]).T
df3=(df.iloc[:,0:1]).T
final=(concat([df3,df2]).T)


# remove identical line
bases=['A','C','G','T']


for i in bases:

    final=final[final !=i]  #!=i
    final=final.dropna(how='all')
    final=final.fillna(i)


reversed=final.T


#save file with output name for fasta
with open(output_file,'w') as output:
    reversed.to_csv(output, sep='\t')

#save total file for plotting
with open(output2_file,'w') as output2:
    final.to_csv(output2, sep='\t')


