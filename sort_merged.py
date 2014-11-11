#for replacement of a given value with NaN
#http://stackoverflow.com/questions/18172851/deleting-dataframe-row-in-pandas-based-on-column-value

# to remove any symbol and replace it with nan
#http://stackoverflow.com/questions/875968/how-to-remove-symbols-from-a-string-with-python


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

parser.add_argument('-o', '--output', help="sorted snp table for fasta conversion")
parser.add_argument('-s', '--snp_table', help="snp table to sort")
parser.add_argument('-t', '--total', help="inverted table to look at")

args = parser.parse_args()
output_file = args.output
input_file = args.snp_table
output2_file = args.total

#read in file as dataframe
df = read_csv(input_file, index_col=[0,1], header=[0,1], dtype=unicode)
df=df.drop('syn/nsyn/intergenic', axis=1 , level=1) 

#replaces lines with "No Hits" with NaN and removes lines with NaN in qbase columns
no_hit= df.mask(df=='No Hit')
removed=no_hit.dropna()


#replaces lines with indel

indel=removed.mask(removed=='indel')
indel2=indel.dropna()



#slash=indel2.replace('/.*', 'NaN')
#slash2=slash.dropna()


# remove identical line
bases=['A','C','G','T']


for i in bases:

    indel2=indel2[indel2 !=i]  #!=i
    indel2=indel2.dropna(how='all')
    indel2=indel2.fillna(i)
#i


#creates dataframe with rows that have / in it
count_qbase=indel2.columns.values
qindexes=[]
for i, v in enumerate(count_qbase):
    if 'qbase:' in v[1]:
        qindexes.append(i)

empty=DataFrame()

for x in qindexes:
   
    empty=empty.append((indel2[(indel2.iloc[:,x]).str.contains('/')]))

t=1+len(qindexes)
# remove rows that are in empty
final=indel2.drop(empty.index[:])
final2 = (final.reset_index(drop=True)).T #drops indexes
final3=final2.iloc[0:t,:]
final3=final3.reset_index() #puts back indexes into dataframe
final3=final3.drop([final3.columns[0]], axis=1) #removes galaxy column




#save file with output name for fasta
with open(output_file,'w') as output:
    final3.to_csv(output, sep='\t')

#save total file for sorting and plotting
with open(output2_file,'w') as output2:
    final.to_csv(output2, sep='\t')


