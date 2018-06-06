import sys
import tempfile
import pandas as pd
import os

fileName = sys.argv[1]
AdpSig = 0
KmerSig = 0
AdapterFile = open("tmp1","a")
Kmer_count = 0

def trigStart_adapter(line):
        global AdpSig
        if "Adapter Content" in line:
                AdpSig = 1
                adapterJudgement = line.rstrip().split("\t")[1]
                print("adapterJudgement: %s" %(adapterJudgement))
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
                kmerJudgement = line.rstrip().split("\t")[1]
                print("kmerJudgement: %s" %(kmerJudgement))
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
                if AdpSig and "Adapter Content" not in line:
                        AdapterFile.write(line)
                if KmerSig and "Kmer Content" not in line:
                        Kmer_count += 1
        AdapterFile.close()

if len(open("tmp1").readlines()) == 1:
        print("no adapter left")
else:
        df = pd.read_table("tmp1",sep="\t",index_col=0)
        df["sum"] = df.sum(axis=1)
        print(max(df["sum"]))
os.remove("tmp1")
print("Kmer_count: %s" %(Kmer_count))                                     
