import sys
import pandas as pd
from pandas import DataFrame,Series
import numpy as np

f1=open(sys.argv[1]).readlines()     #strain_list
f2=sys.argv[2]     #loci_file_list
strainList = [i.rstrip() for i in f1]
lociFrame = pd.read_table(f2,names=["loci","ref"],index_col="loci")
refDict = dict(zip(lociFrame.index,lociFrame["ref"]))
locationList = lociFrame.index
strainFas = []
outputFas = open("phylogeny.fas","a")

for strain in strainList:
        snpFile = open(strain+".fas").readlines()[1]
        snpList = list(snpFile.rstrip().upper())
        strainFas.append(snpList)

frame = pd.DataFrame(strainFas,index=strainList,columns=locationList)
locationListModified = []

for loci in locationList:
        outcome = frame[loci].unique()
        if "N" in outcome and len(outcome) == 2: 
                pass
        elif len(outcome) == 1:
                pass
        else:
                locationListModified.append(loci)

newFrame = frame[locationListModified]
lociFrame = lociFrame.ix[locationListModified]

for strain in strainList:
        outputFas.write(">%s\n" %(strain))
        temp = newFrame.ix[strain]
        outputFas.write("%s\n" %("".join(temp)))

def processSnpFile(strain,locationList,path=None):
        if path == None:
                snpFrame = pd.read_table(strain+".snp",names=["loci","ref","mut"],index_col="loci")
        else:
                snpFrame = pd.read_table(path+"/"+strain+".snp",names=["loci","ref","mut"],index_col="loci")
        fasList = []
        snpDict = dict(zip(snpFrame.index,snpFrame["mut"]))
        for loci in locationList:
                if snpDict.get(loci):
                        fasList.append(snpDict[loci])
                else:
                        fasList.append(refDict[loci])
        return fasList

outputFas.write(">%s\n" %("canettii"))
outputFas.write("%s\n" %("".join(processSnpFile("canettii",locationListModified,path="/home/zty"))))

lociFrame.insert(0,"ref_copy",lociFrame.index)
lociFrame.to_csv("lociFileModified.list",sep="\t",header=False)
