##
## run as: python ace.py --file Analysis/closed_otupicking/otu_table.csv -rt 4
## returns: a two-column csv file with sample ID and estimated ACE value
##
import argparse
import skbio as sk
import numpy as np
import pandas as pd
import skbio.diversity.alpha


parser = argparse.ArgumentParser(description='Estimate ACE.')
parser.add_argument('-f','--file', metavar='FILE',
                   help='File with OTU counts (external OTU table)')
parser.add_argument('-rt','--rare_threshold', metavar='THRESHOLD', type=int,
                   help='Threshold to consider an OTU rare (N. of samples)')

args = parser.parse_args()

print("OTU Table")
print(args.file)
print("Threshold")
print(args.rare_threshold)

#read tab separated otu table
otu= pd.read_csv(args.file,sep="\t",skiprows=1)
otu=otu.rename(columns = {'#OTU ID':'otu'}) #change name of first column

#skbio.diversity.alpha.ace(counts, rare_threshold=args.rare_threshold)

ace = open('ace.csv','w')

for cols in otu.drop(['otu','taxonomy'],axis=1):
        counts= np.array(otu[cols]).astype(np.int64)
        print(skbio.diversity.alpha.ace(counts, rare_threshold=args.rare_threshold))
        ace.write('%s,%s\n' % (cols, skbio.diversity.alpha.ace(counts, rare_threshold=args.rare_threshold)))

ace.close()
print("ACE values written to file ace.csv")
