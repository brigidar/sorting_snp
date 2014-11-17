#for replacement of a given value with NaN
#http://stackoverflow.com/questions/18172851/deleting-dataframe-row-in-pandas-based-on-column-value

# to remove any symbol and replace it with nan
#http://stackoverflow.com/questions/875968/how-to-remove-symbols-from-a-string-with-python

# for isin information
#http://pandas.pydata.org/pandas-docs/stable/indexing.html

#!/usr/bin/env python


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
df=df.drop('syn/nsyn/intergenic', axis=1)

#replaces lines with "No Hits" with NaN and removes lines with NaN in qbase columns
no_hit= df.mask(df=='No Hit')
removed=no_hit.dropna()

# only columns with qbase and refbase in table
count_qbase=list(removed.columns.values)
qindexes=[]
for i, v in enumerate(count_qbase):
    if 'qbase:' in v:
        qindexes.append(i)
df2=(removed.iloc[:,qindexes]).T
df3=(removed.iloc[:,0:1]).T
df4=(concat([df3,df2]).T)

#replaces lines with indel

indel=df4.mask(removed=='indel')
indel2=indel.dropna()


# remove identical line
bases=['A','C','G','T']


for i in bases:

    indel2=indel2[indel2 !=i]  #!=i
    indel2=indel2.dropna(how='all')
    indel2=indel2.fillna(i)



#creates dataframe with rows that have / in it


empty=DataFrame()
t=1+len(qindexes)
for x in qindexes:
   
    empty=empty.append((indel2[(indel2.iloc[:,x]).str.contains('/')]))


# remove rows that are in empty
final=(indel2.drop(empty.index[:]))
final2 = (final.reset_index(drop=True)).T #drops indexes


#description rows back in overview table look at isin options with index it checks if the index is included in the other index and only keeps the one that are
i=t+1
rest=df.iloc[:,i:]
sel=rest[rest.index.isin(final.index)]
sel=sel.T
final=final.T
final3=(concat([final, sel])).T


#save file with output name for fasta -o option
with open(output_file,'w') as output:
    final2.to_csv(output, sep='\t')

#save total file for plotting -t option
with open(output2_file,'w') as output2:
    final3.to_csv(output2, sep='\t')


