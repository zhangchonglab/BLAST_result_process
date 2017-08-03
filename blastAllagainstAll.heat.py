# this script is used to reformat the strucutre of outfmt 6 of blast all against all for het map drawing
# input allagaintall outfmt 6

import os
import sys
blastFile=sys.argv[1]
queryLst=sys.argv[2]

QueryLst=[]
AllvsAll_dic={}
f=open(blastFile,'r')
k=open(queryLst,'r')
for line in k:
    gene=line.rstrip()
    QueryLst.append(gene)
    AllvsAll_dic[gene]={}
k.close()

# the basic dic to store all vs all data of identity
for gene in AllvsAll_dic:
    for hit in QueryLst:
        AllvsAll_dic[gene][hit]=0.0

# read the data in
for line in f:
    row=line.rstrip().split('\t')
    query=row[0]
    hit=row[1]
    identity=row[2]
    AllvsAll_dic[query][hit]=identity
f.close()

# write the data to a flatfile
os.system('cat /dev/null > %s.reformat.txt'%(blastFile))
g=open('%s.reformat.txt'%(blastFile),'r+')
# the first line
writtenline='data\t'
for query in QueryLst:
    writtenline+=query+'\t'
g.write(writtenline[:-1]+'\n')

for query in QueryLst:
    writtenline='%s\t'%(query)
    for hit in QueryLst:
        writtenline+=AllvsAll_dic[query][hit]+'\t'
    g.write(writtenline[:-1]+'\n')
g.close()