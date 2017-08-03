# this script is used to extract the desired hits from a blast result based on several thresholds:
# identity, coverage for query, coverage for hit;
# utility: python blast_filter.py query db blast_result
# query: query fasta file used in the abovementioned BLAST search
# db: databse fasta file used in the abovementioned BLAST search
# blast_result output file of BLAST+ package (BLASTP, BLASTN, etc) in outfmt6 format

# The program will ask the user to type in the thresholds (identity, coverage for query, coverage for hit) used for hit extraction. 
# After run the script, you can see tips for typing in, just follow it, please.
# e.g. "Whether you would like to use identity as threshold (yes or no):"
# "Type in the cutoff for identity (0.0-1.0 real number):"

# first of all get the query and database fastafile, respectively
import os
import sys

query_file=sys.argv[1]
database_file=sys.argv[2]
blast_result_file=sys.argv[3]

# first of all, ask the user to define the threshold that they want to use:
# identity
flag='no'
identity_thre=False
identity=0.0
flag=raw_input("Whether you would like to use identity as threshold (yes or no):   ")
if flag=='yes':
    identity_thre=True
    identity=float(raw_input("Type in the cutoff for identity (0.0-1.0 real number):   "))
print 'Identity: %s'%str(identity)

# coverage for query
flag='no'
queryCoverage_thre=False
queryCoverage=0.0
flag=raw_input("Whether you would like to use query coverage as threshold (yes or no):   ")
if flag=='yes':
    queryCoverage_thre=True
    queryCoverage=float(raw_input("Type in the cutoff for query coverage (0.0-1.0 real number):   "))
print 'queryCoverage: %s'%str(queryCoverage)

# coverage for hit
flag='no'
hitCoverage_thre=False
hitCoverage=0.0
flag=raw_input("Whether you would like to use hit coverage as threshold (yes or no):   ")
if flag=='yes':
    hitCoverage_thre=True
    hitCoverage=float(raw_input("Type in the cutoff for query coverage (0.0-1.0 real number):   "))
print 'hitCoverage: %s'%str(hitCoverage)

# now construct the query and hit length dic with data structure like this
# {queryID:float,queryID:float}
# the same as hit
# query_length_dic and hit_length_dic, respectively
query_length_dic={}
proteinID='' # to store the protein ID now
f=open(query_file,'r')
for line in f:
    if line[0]=='>':
        proteinID=line[1:(-1)].split(' ')[0]
        query_length_dic[proteinID]=0.0
    else:
        query_length_dic[proteinID]+=float(len(line.rstrip()))
f.close()        

hit_length_dic={}
proteinID='' # to store the protein ID now
f=open(database_file,'r')
for line in f:
    if line[0]=='>':
        proteinID=line[1:(-1)].split(' ')[0]
        hit_length_dic[proteinID]=0.0
    else:
        hit_length_dic[proteinID]+=float(len(line.rstrip()))
f.close()        

# now it is time for filtering
f=open(blast_result_file,'r')
os.system('cat /dev/null > %s.trim.iden%s.queCo%s.hitCo%s'%(blast_result_file,str(int(identity*100)),str(int(queryCoverage*100)),str(int(hitCoverage*100))))
g=open('%s.trim.iden%s.queCo%s.hitCo%s'%(blast_result_file,str(int(identity*100)),str(int(queryCoverage*100)),str(int(hitCoverage*100))),'r+')
keys=['queryId', 'subjectId', 'percIdentity', 'alnLength', 'mismatchCount', 'gapOpenCount', 'queryStart', 'queryEnd', 'subjectStart', 'subjectEnd', 'eVal', 'bitScore']
blast_result=[]
for line in f:
    row=dict.fromkeys(keys)
    elements=line.rstrip().split('\t')
    for i, elem in enumerate(elements):
        row[keys[i]]=elem
    flag=True  # judge whether this result pass or not
    if identity_thre:
        if float(row['percIdentity'])/100.0<identity:
            flag=False
    if queryCoverage_thre:
        result_queryCoverage=float(row['alnLength'])/float(query_length_dic[row['queryId']])
        if result_queryCoverage<queryCoverage:
            flag=False
    if hitCoverage_thre:
        result_hitCoverage=float(row['alnLength'])/float(hit_length_dic[row['subjectId']])
        if result_hitCoverage<hitCoverage:
            flag=False
    if flag:
        g.write(line)
f.close()
g.close()
