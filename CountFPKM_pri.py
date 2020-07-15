#Author: SunJing
#Date: 2020.07.15
#Description: Count pri-miRNA FPKM from htseq-count result.
import sys,re
usage="""
********************************************************************
Usage:
    python CountFPKM_pri.py 
    
input files:
    pri-mirna.gff: 325 pri-miRNAs annotated in Araport11
    total mapped reads file: tab seperated; including headline; total mapped reads were put in the last column
    pri-miRNA counts file: tab seperated; merge the htseq-count result of each lib together; sample order must be the same as total mapped reads file
********************************************************************
"""

if len(sys.argv) > 1:
    print usage
    sys.exit(0)

def main():
    idfile = open('pri-mirna.gff','r')
    a = idfile.readline()
    mrna = {}
    name = {}
    for line in idfile:
        if line[0] == '#':
            continue
        rec = line.strip('\n').split('\t')
        ID = rec[-1].split(';')[0].strip('ID=')
        mrna[ID] = int(rec[4])-int(rec[3])+1
        name[ID] = rec[-1].split(';')[3].strip('Name=')
    idfile.close()
    print mrna
    idfile = open('data_info_root.txt','r')
    a = idfile.readline()
    readsnum = {}
    i = 0
    for line in idfile:
        rec = line.split('\t')
        i += 1
        readsnum[i] = float(rec[-1])
    idfile.close()
    
    idfile = open('root-pri_counts.txt','r')
    a = idfile.readline()
    out = open('root-pri.fpkm','w')
    out.write(a)
    for line in idfile:
        rec = line.split('\t')
        if rec[0] not in mrna:
            continue
        out.write(rec[0])
        for i in range(1,len(rec)):
            FPKM = float(rec[i])/mrna[rec[0]]/readsnum[i]*1000000000
            out.write('\t'+str(FPKM))
        out.write('\t'+name[rec[0]]+'\n')
    idfile.close()
    out.close()


if __name__ == "__main__":
    main()

