import sys
import pandas as pd
from pandas import DataFrame,Series

f1=open(sys.argv[1]).readlines()     #strain_list
f2=sys.argv[2]     #loci_file_list
strainList = [i.rstrip() for i in f1]
lociFrame = pd.read_table(f2,names=["loci","ref"],index_col="loci")
refDict = dict(zip(lociFrame.index,lociFrame["ref"]))
locationList = lociFrame.index
strainFas = []
outputFas = open("phylogeny.fas","a")

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

for strain in strainList:
        strainFas.append(processSnpFile(strain,locationList))

frame = pd.DataFrame(strainFas,index=strainList,columns=locationList)
locationListModified = []

for loci in locationList:
        outcome = frame[loci].unique()
        if len(outcome) == 1:
                print("%s is invariant" %(loci))
        else:
                locationListModified.append(loci)

newFrame = frame[locationListModified]
lociFrame = lociFrame.ix[locationListModified]

for strain in strainList:
        outputFas.write(">%s\n" %(strain))
        temp = newFrame.ix[strain]
        outputFas.write("%s\n" %("".join(temp)))

outputFas.write(">%s\n" %("canettii"))
outputFas.write("%s\n" %("".join(processSnpFile("canettii",locationListModified,path="/home/zty"))))

lociFrame.insert(0,"refCopy",lociFrame.index)
lociFrame.to_csv("lociFileModified.list",sep="\t",header=False)
