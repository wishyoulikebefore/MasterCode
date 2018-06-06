import sys
import tempfile
import pandas as pd
import os

fileName = sys.argv[1]
AdpSig = 0
KmerSig = 0
AdapterFile = open("tmp1","a")
KmerFile = open("tmp2","a")

def trigStart_adapter(line):
        global AdpSig
        if "Adapter Content" in line:
                AdpSig = 1
        return

def trigEnd_adapter(line):
        global AdpSig
        if "END_MODULE" in line and AdpSig:
                AdpSig = 0
        return

def trigStart_kmer(line):
        global KmerSig
        if "Kmer Content" in line:
                KmerSig = 1
        return

def trigEnd_kmer(line):
        global KmerSig
        if "END_MODULE" in line and KmerSig:
                KmerSig = 0
        return

with open(fileName) as f:
        for line in f.readlines():
                trigStart_adapter(line)
                trigEnd_adapter(line)
                trigStart_kmer(line)
                trigEnd_kmer(line)
                if AdpSig:
                        AdapterFile.write(line)
                if KmerSig:
                        KmerFile.write(line)
        AdapterFile.close()
        KmerFile.close()

fastqc_judgement = open("tmp1").readlines()[0].rstrip().split("\t")[1]
print(fastqc_judgement)
if len(open("tmp1").readlines()) == 3:
        print("no adapter left")
else:
        dflist = []
        columnsName = open("tmp1").readlines()[1].rstrip().split("\t")[1:]
        indexName = []
        for line in open("tmp1").readlines()[2:-1]:
                lineOfList = line.rstrip().split("\t")
                indexName.append(lineOfList[0])
                dflist.append([float(i) for i in lineOfList[1:]])
        df = pd.DataFrame(dflist,columns=columnsName,index=indexName)
        df["sum"] = df.sum(axis=1)
        print(max(df["sum"]))
os.remove("tmp1")

kmer = len(open("tmp2").readlines())
print(kmer-2)
os.remove("tmp2")
