import sys
from pandas import DataFrame,Series
import pandas as pd
import numpy as np
import multiprocessing

f1 = sys.argv[1]    #markedLineage.list
f2 = sys.argv[2]     #loci_file_list
lociFrame = pd.read_table(f2,names=["loci","ref"],index_col="loci")
lineageFrame = pd.read_table(f1,names=["strain","lineage"])
lineageDict = dict(zip(lineageFrame["strain"],lineageFrame["lineage"]))
lociList = list(lociFrame.index)
strainList = lineageFrame["strain"]
df = DataFrame([],columns=lociList+["lineages"])
df.to_csv("raw_data_MP",index=False)
output = open("raw_data_MP","a")

def processSnpFile(strain):
        print("Start analyse %s" %(strain))
        snpFrame = pd.read_table(strain+".snp",names=["loci","ref","mut"],index_col="loci")
        infoList = []
        snpDict = dict(zip(snpFrame.index,snpFrame["mut"]))
        for loci in lociList:
                if snpDict.get(loci):
                        infoList.append(snpDict[loci])
                else:
                        infoList.append("")
        infoList.append(lineageDict.get(strain))
        print("%s is finished" %(strain))
        return infoList

def mycallback(x):
        output.write(",".join(x)+"\n")
        
if __name__ == "__main__":
        pool = multiprocessing.Pool(processes=24)
        for strain in strainList:
                pool.apply_async(processSnpFile,(strain,),callback=mycallback)
        pool.close()
        pool.join()
