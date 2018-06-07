"""
比较trim前后的fastqc结果
"""

import sys
import os
import pandas as pd

"""
selected file format:
P18006728LU01_combined_R1_fastqc        out_P18006728LU01_combined_R1_fastqc
"""

compareDict = {}
AdpSig = 0
KmerSig = 0
Kmer_count = 0
adapterJudgement = ""
kmerJudgement = ""

def trigStart_adapter(line):
        global AdpSig,adapterJudgement
        if "Adapter Content" in line:
                AdpSig = 1
                adapterJudgement = line.rstrip().split("\t")[1]
        return

def trigEnd_adapter(line):
        global AdpSig
        if "END_MODULE" in line and AdpSig:
                AdpSig = 0
        return

def trigStart_kmer(line):
        global KmerSig,kmerJudgement
        if "Kmer Content" in line:
                KmerSig = 1
                kmerJudgement = line.rstrip().split("\t")[1]
        return

def trigEnd_kmer(line):
        global KmerSig
        if "END_MODULE" in line and KmerSig:
                KmerSig = 0
        return

with open(sys.argv[1]) as f:
        for line in f.readlines():
                lineOfList = line.rstrip().split("\t")
                compareDict[lineOfList[0]] = lineOfList[1]

def process(fileName):
        global Kmer_count,adapterJudgement,kmerJudgement
        AdapterFile = open("tmp","a")
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
        if len(open("tmp").readlines()) == 1:
                adapter_maximum = 0
        else:
                df = pd.read_table("tmp",sep="\t",index_col=0)
                df["sum"] = df.sum(axis=1)
                adapter_maximum = max(df["sum"])
                max_loc = df["sum"].argmax()
        os.remove("tmp")
        this_count = Kmer_count
        Kmer_count = 0
        return adapterJudgement,kmerJudgement,adapter_maximum,this_count,max_loc

for raw,trim in compareDict.items():
        rawFile = os.path.join(raw,"fastqc_data.txt")
        trimmedFile = os.path.join(trim,"fastqc_data.txt")
        raw1,raw2,raw3,raw4,raw5 = process(rawFile)
        trim1,trim2,trim3,trim4,trim5 = process(trimmedFile)
        print("Raw data: %s ; Trimmed data: %s" %(rawFile,trimmedFile))
        print("Adapter: %s --> %s" %(raw1,trim1))
        print("Adapter maximum percentage and its location: %s in %s --> %s in %s" %(raw3,raw5,trim3,trim5))
        print("Kmer: %s --> %s" %(raw2,trim2))
        print("Reserved Kmer count: %s --> %s" %(raw4,trim4))
