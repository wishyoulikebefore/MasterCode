import matplotlib.pyplot as plt
from pandas import DataFrame,Series
import matplotlib as mpl
import numpy as np
import argparse
from collections import defaultdict,Counter
import re
import os

parser = argparse.ArgumentParser("plot drug's proportion")
parser.add_argument("drug",choices={"INH","RIF","EMB","PZA"},help="input target drug")
parser.add_argument("integers",type=int,nargs="+",help="strain number of DR, MDR, pre-XDR, XDR")
args = parser.parse_args()

period = ["DR","MDR","pre-XDR","XDR"]
strainNumOfPeriod = dict(zip(period,args.integers))
mutationDict = defaultdict(list)
offset = [-0.3,-0.1,0.1,0.3]

fig = plt.figure()
ax = plt.gca()
mpl.rcParams["font.size"] = 8
ax.tick_params(top='off',bottom='off',left='on',right='off')
colors = plt.cm.hot(np.linspace(0.9,0.1,4))
maxProportionList = []

def process(drug):
    file = open(args.drug+"_mutation_spectrum")
    for row in file.readlines():
        rowOfList = row.rstrip().split()
        mutationList = re.sub(r'(\[|\]|\')',"",rowOfList[2]).split("+")
        for mutation in mutationList:
            if mutation != "":
                mutationDict[mutation].append(rowOfList[1])
            else:
                pass
    frame = dict2frame(mutationDict)
    frame["sum"] = frame["DR"] + frame["MDR"] + frame["pre-XDR"] + frame["XDR"]
    filterFrame = frame[frame["sum"] > 10].sort_values(by="sum",ascending=False)
    indexList = filterFrame.index
    myplot(filterFrame,indexList)

def dict2frame(mutationDict):
    index = []
    columns = ["DR","MDR","pre-XDR","XDR"]
    calculation = []
    for key in mutationDict.keys():
        counter = Counter(mutationDict.get(key))
        index.append(key)
        outcome = [counter[i] for i in columns]
        calculation.append(outcome)
    frame = DataFrame(calculation,columns=columns,index=index)
    if os.path.exists("%s_proportion" %(args.drug)):
        print("%s_proportion has been updated" %(args.drug))
        os.remove("%s_proportion" %(args.drug))
    frame.to_csv("%s_proportion" %(args.drug))
    return frame

def myplot(filterFrame,indexList):
    for nu in range(len(period)):
        colName = period[nu]
        x_shift = offset[nu]
        barColor = colors[nu]
        frameslice = filterFrame[colName]
        proportion = frameslice/strainNumOfPeriod.get(colName)*100
        xloc = np.arange(len(frameslice))+1+x_shift
        plt.bar(xloc,proportion,0.2,color=barColor,label=colName)
        maxProportionList.append(max(proportion))
    plt.text(0.05,0.95,args.drug,transform = ax.transAxes,fontsize=16,fontweight="bold",va="top")
    plt.ylabel("比例(%)", fontsize=15)
    plt.ylim(0, max(maxProportionList) + 5)
    plt.xlim(0, len(indexList) + 2)
    plt.xticks(range(1, len(indexList) + 1), indexList)
    plt.title("%s耐药突变在不同耐药阶段所占比例" %(args.drug), fontsize=20)
    fig.autofmt_xdate()
    plt.legend()
    plt.show()

process(args.drug)
