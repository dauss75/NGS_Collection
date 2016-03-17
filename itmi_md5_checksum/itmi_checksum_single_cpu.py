#!/usr/bin/python

import sys, glob, os, re, timeit
from subprocess import call, Popen, PIPE

'''
Usage:
python itmi_checksum.py [Top Directory Path] 
runs md5 checksum on bam, vcf, tbi files
last updated: March 15, 2016
'''

def md5checksum(fin):
    tmpFile=reduce(lambda x,y: x.extend(y),fin)
    cmd = "md5sum", "%s" % tmpFile
    tag=Popen(cmd, stdout=PIPE, stdin=PIPE).stdout.read()
    return tag.split(' ')
 
def main():
    f2 = open('summary_single_cpu.txt', 'w')
    f2.write('job starts\n')
    start = timeit.default_timer()
    mainDirPath=sys.argv[1]
    folderList=glob.glob(mainDirPath + "*")
    check={}; checksum={}
    for i in xrange(0, len(folderList)):
        bamFile=glob.glob(folderList[i]+"/Assembly/*.bam")
        bamTag=md5checksum(bamFile)
        bamHashcode=bamTag[0]
        bamName=re.match(r".*/(Assembly.*.bam)",bamTag[2]).group(1)
        check[bamName]=bamHashcode          

        snpFile=glob.glob(folderList[i]+"/Variations/*.vcf.gz")
        snpTag=md5checksum(snpFile)
        snpHashcode=snpTag[0]
        snpName=re.match(r".*/(Variations.*.vcf.gz)",snpTag[2]).group(1)
        check[snpName]=snpHashcode        

        tbiFile=glob.glob(folderList[i]+"/Variations/*.tbi")
        tbitag=md5checksum(tbiFile)
        tbiHashcode=tbitag[0]
        tbiName=re.match(r".*/(Variations.*.tbi)",tbitag[2]).group(1)
        check[tbiName]=tbiHashcode
        
#       open md5checksum file
        checksumFile=glob.glob(folderList[i]+"/*txt*")       
        checksumFile=reduce(lambda x,y: x.extend(y),checksumFile)
        
        f1=open(checksumFile, 'r')
        for line in f1:
            listHash=line.strip().split(' ')[0]
            listFile=line.strip().split(' ')[2]
            checksum[listFile]=listHash
        
#       check common items and hash codes
        shared_items = set(check.items()) & set(checksum.items())
        if len(shared_items) != 3:
            f2.write(folderList[i]+ " fail\n")
        else:
            f2.write(folderList[i]+ " pass\n")
            
        check.clear(); checksum.clear()     
    stop = timeit.default_timer()
    runTime=stop - start
    f2.write('job done!\n')
    f2.write("elasped time: %.2fsec" % runTime)
    f1.close(); f2.close()
if __name__ == '__main__':
    main()
