#Author: SunJing
#Date: 2020.07.15
#Description: Count FPKM from htseq-count result.
import sys,re
usage="""
********************************************************************
Usage:
    python CountFPKM.py 
    
input files:
    Araport11.gff
    total mapped reads file: tab seperated; including headline; total mapped reads were put in the last column
    read counts file: tab seperated; merge the htseq-count result of each lib together; sample order must be the same as total mapped reads file
********************************************************************
"""

if len(sys.argv) > 1:
    print usage
    sys.exit(0)
def main():
    idfile = open('/disk2/jsun/reference/TAIR/Araport11.gff','r')
    a = idfile.readline()
    mrna = {}
    ID = ''
    mrna_ID = ''
    gene_ID = ''
    for line in idfile:
        if line[0] == '#':
            continue
        rec = line.strip('\n').split('\t')
        #if rec[2] == 'CDS' or rec[2] == 'five_prime_UTR' or rec[2] == 'three_prime_UTR'
        if rec[2] == 'exon':
            if rec[-1].split(';')[1].replace('Parent=','').split('.')[0] != gene_ID:
                #print gene_ID,mrna_ID
                if gene_ID:
                    cdslist.append(cds_len)
                    mrna[gene_ID] = max(cdslist)
                cdslist = []
                gene_ID = rec[-1].split(';')[1].replace('Parent=','').split('.')[0]
                mrna_ID = rec[-1].split(';')[1].replace('Parent=','')
                cds_len = 0
            else:
                if rec[-1].split(';')[1].replace('Parent=','') != mrna_ID:
                    cdslist.append(cds_len)
                    cds_len = 0
                    mrna_ID = rec[-1].split(';')[1].replace('Parent=','')
            cds_len += int(rec[4]) - int(rec[3]) + 1
    cdslist.append(cds_len)
    mrna[gene_ID] = max(cdslist)
    idfile.close()

    idfile = open('data_info_root.txt','r')
    a = idfile.readline()
    readsnum = {}
    i = 0
    for line in idfile:
        rec = line.split('\t')
        i += 1
        readsnum[i] = float(rec[-1])
    idfile.close()
    
    idfile = open('root_counts.txt','r')
    a = idfile.readline()
    out = open('root.fpkm','w')
    out.write(a)
    for line in idfile:
        rec = line.split('\t')
        if rec[0] not in mrna:
            continue
        out.write(rec[0])
        for i in range(1,len(rec)):
            FPKM = float(rec[i])/mrna[rec[0]]/readsnum[i]*1000000000
            out.write('\t'+str(FPKM))
        out.write('\n')
    idfile.close()
    out.close()


if __name__ == "__main__":
    main()

