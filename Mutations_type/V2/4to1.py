"""
最大的问题是autofmt_xdate()会关闭bottom以外subplots的X轴
"""

import matplotlib.pyplot as plt
from pandas import DataFrame,Series
import matplotlib as mpl
import numpy as np
import argparse
from collections import defaultdict,Counter
import re
import os
from pprint import pprint

parser = argparse.ArgumentParser("plot drug's proportion")
parser.add_argument("drugs",nargs="+",choices={"INH","RIF","EMB","PZA"},help="input target drug")
args = parser.parse_args()

period = ["DR","MDR","pre-XDR","XDR"]
offset = [-0.3,-0.1,0.1,0.3]
colors = plt.cm.hot(np.linspace(0.9,0.1,4))

def process(drug):
    file = open(drug+"_mutation_spectrum")
    for row in file.readlines():
        rowOfList = row.rstrip().split()
        mutationList = re.sub(r'(\[|\]|\')',"",rowOfList[2]).split("+")
        for mutation in mutationList:
            if mutation != "":
                mutationDict[mutation].append(rowOfList[1])
            else:
                pass
    frame = dict2frame(mutationDict,drug)
    frame["sum"] = frame["DR"] + frame["MDR"] + frame["pre-XDR"] + frame["XDR"]
    filterFrame = frame[frame["sum"] > 10].sort_values(by="sum",ascending=False)
    indexList = filterFrame.index
    file.close()
    return filterFrame,indexList

def dict2frame(mutationDict,drug):
    index = []
    calculation = []
    for key in mutationDict.keys():
        counter = Counter(mutationDict.get(key))
        index.append(key)
        outcome = [counter[i] for i in period]
        calculation.append(outcome)
    frame = DataFrame(calculation,columns=period,index=index)
    if os.path.exists("%s_proportion" %(drug)):
        print("%s_proportion has been updated" %(drug))
        os.remove("%s_proportion" %(drug))
    frame.to_csv("%s_proportion" %(drug))
    return frame

def myplot(filterFrame,indexList,ax):
    for nu in range(len(period)):
        colName = period[nu]
        x_shift = offset[nu]
        barColor = colors[nu]
        frameslice = filterFrame[colName]
        proportion = frameslice/strainNumEachPeriod.get(colName)*100
        xloc = np.arange(len(frameslice))+1+x_shift
        ax.bar(xloc,proportion,0.2,color=barColor,label=colName)
        maxProportionList.append(max(proportion))
    ax.text(0.1,1.1,drug,transform = ax.transAxes,fontsize=16,fontweight="bold",va="top")
#    plt.ylim(0, max(maxProportionList) + 5)
#    plt.xlim(0, len(indexList) + 2)
    ax.set_xticks(range(1,len(indexList)+1))
    ax.set_xticklabels(indexList)
    if i == 0:
        plt.legend()

fig = plt.figure(figsize=(16,10))
ax = plt.gca()
fig.suptitle("IREP_plot",fontsize=16,fontstyle="italic",fontweight="bold")
fig.autofmt_xdate()
plt.axis("off")
plt.text(-0.1,0.6,"比例(%)",transform = ax.transAxes,fontsize=16,rotation=90)

for i,drug in enumerate(args.drugs):
    ax = fig.add_subplot(2,2,i+1)
    mutationDict = defaultdict(list)
    filterFrame,indexList = process(drug)
    maxProportionList = []
    strainNumEachPeriod = {}
    if drug in ["INH","RIF"]:
        strainNumEachPeriod["DR"] = len(open("DR_"+drug+".list").readlines())
        strainNumEachPeriod["MDR"] = len(open("MDR.list").readlines())
        strainNumEachPeriod["pre-XDR"] = len(open("pre-XDR.list").readlines())
        strainNumEachPeriod["XDR"] = len(open("XDR.list").readlines())
    else:
        for periodName in period:
            strainNumEachPeriod[periodName] = len(open(periodName+"_"+drug+".list").readlines())
    myplot(filterFrame,indexList,ax)

plt.show()

