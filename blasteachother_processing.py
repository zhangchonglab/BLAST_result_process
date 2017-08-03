# this sciprt is used to restructure the outfmt 6 blast each other result in order to draw a heatmap 
# input data strucutre:
# outfmt 6 blast result
# output data set:
# NAME  PROT1 PROT2 ...
# PROT1 xxxx  xxxx  ...
# PROT2 XXXX  XXXX  ...
# ...

import os
import sys
blastresult=sys.argv[1]
Termspecified=int(sys.argv[2])-1

# the second term refer to the term you want to add as variable here,
# 3 for identity and 4 for queryCoverage

f=open(blastresult,'r')
protein_lst=[]
protein_term_dic={}

for line in f:
    row=line.rstrip().split('\t')
    protein=row[0]
    if protein not in protein_lst:
        protein_lst.append(protein)
        protein_term_dic[protein]={}
f.close()

protein_lst.sort()

for protein in protein_lst:
    for hit in protein_lst:
        protein_term_dic[protein][hit]=0

f=open(blastresult,'r')
for line in f:
    row=line.rstrip().split('\t')
    protein=row[0]
    hit=row[1]
    term=row[Termspecified]
    # for the order of the hit is ranked by evalue
    if protein_term_dic[protein][hit]==0:
        protein_term_dic[protein][hit]=float(term)
f.close()

if Termspecified==3:
    for protein in protein_lst:
        max=0
        for hit in protein_lst:
            if protein_term_dic[protein][hit] > max:
                max=protein_term_dic[protein][hit]
        for hit in protein_lst:
            protein_term_dic[protein][hit]=protein_term_dic[protein][hit]/max*100
else:
    print 'not query coverage'

print 'Constructed'

os.system('cat /dev/null > %s.heatmap.%d.txt'%(blastresult,Termspecified))
g=open('%s.heatmap.%d.txt'%(blastresult,Termspecified),'r+')
written_line='NAME\t'
for protein in protein_lst:
    written_line+=protein+'\t'
g.write(written_line[:-1]+'\n')

for protein in protein_lst:
    written_line=protein+'\t'
    for hit in protein_lst:
        written_line+=str(protein_term_dic[protein][hit])+'\t'
    g.write(written_line[:-1]+'\n')
g.close()


