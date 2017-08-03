# this script utilized a blastn/p outfmt6 file and database file to extract the exactly mapped seq along the query according to the BLAST search result.
# Several constant nucleotides (amino acids) can be also extracted, the number of which is specified by upstream_number and downstream_number
# utility: python MappingExtract.py blast_result db upstream_number downstream_number prefix
# blast_result output file of BLAST+ package (BLASTP, BLASTN, etc) in outfmt6 format
# db: databse fasta file used in the abovementioned BLAST search
# upstream_number: number of nucl or aa upstream of query, integer
# downstream_number: number of nucl or aa downstream of query, integer
# prefix: prefix to name the output files

import os
import sys
from Bio.Seq import Seq

blastResult=sys.argv[1]
fastaFile=sys.argv[2]
Upstream=int(sys.argv[3])
Downstream=int(sys.argv[4])
prefix=sys.argv[5]

# construct the list for the fasta
FastaSeq=''
f=open(fastaFile,'r')
for line in f:
    if line[0]!='>':
        FastaSeq+=line.rstrip()
f.close()
print 'fasta file successfully processed'

# processing the blast result file
f=open(blastResult,'r')
os.system('cat /dev/null > %s.fasta'%(prefix))
g=open('%s.fasta'%(prefix),'r+')
for line in f:
    row=line.rstrip().split('\t')
    queryID=row[0]
    g.write('>%s\n'%(queryID))
    Hit1st=int(row[8])
    Hit2nd=int(row[9])
    if Hit1st<Hit2nd:
        wantedStart=Hit1st-Upstream-1
        wantedEnd=Hit2nd+Downstream
        mappingSeq=FastaSeq[wantedStart:wantedEnd]
        g.write(mappingSeq+'\n')
    else:
        wantedStart=Hit1st+Upstream
        wantedEnd=Hit2nd-Downstream-1
        mappingSeq=Seq(FastaSeq[wantedEnd:wantedStart]).reverse_complement()
        g.write(str(mappingSeq[:])+'\n')
f.close()   
g.close()
