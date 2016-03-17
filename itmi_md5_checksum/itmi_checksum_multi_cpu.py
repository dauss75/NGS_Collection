#!/usr/bin/python

import sys, glob, os, re, timeit
from subprocess import call, Popen, PIPE
from multiprocessing import Pool, cpu_count
'''
Usage:
python itmi_checksum_multi_cpu.py [Top Directory Path] 
runs md5 checksum on bam, vcf, tbi files
last updated: March 15, 2016
'''

check={}; checksum={}
N_CPU=8
N_ITEM=3

def md5checksum(fin):
    tmpFile=reduce(lambda x,y: x.extend(y),fin)
    cmd = "md5sum", "%s" % tmpFile
    tag=Popen(cmd, stdout=PIPE, stdin=PIPE).stdout.read()
    return tag.split(' ')

def process_files(name):
    check={}; checksum={}
    bamFile=glob.glob(name+'/Assembly/*.bam')
    bamTag=md5checksum(bamFile)
    bamHashcode=bamTag[0]
    bamName=re.match(r".*/(Assembly.*.bam)",bamTag[2]).group(1)
    check[bamName]=bamHashcode
    
    snpFile=glob.glob(name+'/Variations/*.vcf.gz')
    snpTag=md5checksum(snpFile)
    snpHashcode=snpTag[0]
    snpName=re.match(r".*/(Variations.*.vcf.gz)",snpTag[2]).group(1)
    check[snpName]=snpHashcode

    tbiFile=glob.glob(name+'/Variations/*.tbi')
    tbitag=md5checksum(tbiFile)
    tbiHashcode=tbitag[0]
    tbiName=re.match(r".*/(Variations.*.tbi)",tbitag[2]).group(1)
    check[tbiName]=tbiHashcode

#   open md5checksum file
    checksumFile=glob.glob(name+'/*txt*')
    checksumFile=reduce(lambda x,y: x.extend(y),checksumFile)

    f1=open(checksumFile, 'r')
    for line in f1:
        listHash=line.strip().split(' ')[0]
        listFile=line.strip().split(' ')[2]
        checksum[listFile]=listHash

#   check common items and hash codes
    shared_items = set(check.items()) & set(checksum.items())
    check.clear(); checksum.clear(); f1.close() 
    if len(shared_items) != N_ITEM:
        status=name+" fail\n"
        return status
    else:
        status=name+" pass\n"
        return status 
def main():
    f2 = open('summary_multiprocess.txt', 'w')
    f2.write('job starts\n')
    start = timeit.default_timer()
    mainDirPath=sys.argv[1]
    folderList=glob.glob(mainDirPath + "*")
#    ncpu=cpu_count()
    ncpu=N_CPU
    pool = Pool(processes=ncpu)
    result= pool.map(process_files, folderList)
    pool.close()
    pool.join()
    f2.write('parallelization done!\n')
    for i in range(len(result)):
        f2.write(result[i])
            
    stop = timeit.default_timer()
    runTime=stop - start
    f2.write('number of cpu: %d\n'  % ncpu)
    f2.write("elasped time: %.2fsec" % runTime)
    f2.close()
if __name__ == '__main__':
    main()
