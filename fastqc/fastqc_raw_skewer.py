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

for sample in open(sys.argv[1]).readlines():
        sample = sample.rstrip()
        raw1 = os.path.join(sample+"_combined_R1_fastqc","fastqc_data.txt")
        raw2 = os.path.join(sample+"_combined_R2_fastqc","fastqc_data.txt")
        trim1 = os.path.join(sample+"_combined_R1_trimmed_fastqc","fastqc_data.txt")
        trim2 = os.path.join(sample+"_combined_R2_trimmed_fastqc","fastqc_data.txt")
        raw1_a,raw1_b,raw1_c,raw1_d,raw1_e = process(raw1)
        raw2_a,raw2_b,raw2_c,raw2_d,raw2_e = process(raw2)
        trim1_a,trim1_b,trim1_c,trim1_d,trim1_e = process(trim1)
        trim2_a,trim2_b,trim2_c,trim2_d,trim2_e = process(trim2)
        print("%s rawdata vs skewer trimmed data" %(sample))
        print("R1")
        print("Adapter: %s --> %s" %(raw1_a,trim1_a))
        print("Adapter maximum percentage and its location: %s in %s --> %s in %s" %(raw1_c,raw1_e,trim1_c,trim1_e))
        print("Kmer: %s --> %s" %(raw1_b,trim1_b))
        print("Reserved Kmer count: %s --> %s" %(raw1_d,trim1_d))
        print("R2")
        print("Adapter: %s --> %s" %(raw2_a,trim2_a))
        print("Adapter maximum percentage and its location: %s in %s --> %s in %s" %(raw2_c,raw2_e,trim2_c,trim2_e))
        print("Kmer: %s --> %s" %(raw2_b,trim2_b))
        print("Reserved Kmer count: %s --> %s" %(raw2_d,trim2_d))
        print()

