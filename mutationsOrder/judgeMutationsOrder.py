import sys

dbFile = open("out.db").readlines()
mutation1Loc = str(sys.argv[1])
mutation1 = sys.argv[2]
mutation2Loc = str(sys.argv[3])
mutation2 = sys.argv[4]
outNode = open("out.node").readlines()

lineNum = len(dbFile)
nodeNum = int((lineNum+1)/6)
mut1 = str(mutation1Loc)+mutation1
mut2 = str(mutation2Loc)+mutation2
mut1EmergeNode = {}
mut2EmergeNode = {}
outNodeDict = {}
mut1InNode = open(mut1+"_node","a")
mut1InTermianl = open(mut1+"_terminalStrain","a")
mut2InNode = open(mut2+"_node","a")
mut2InTermianl = open(mut2+"_terminalStrain","a")
output1 = open("1before2","a")
output2 = open("2before1","a")

for line in outNode:
        lineOfList = line.rstrip().split()
        nodeName = lineOfList[1].replace(">","")
        outNodeDict[nodeName] = lineOfList[0]

def detect_mutation(loc,mutation,nodeNum,rank=1):
        nodeName = dbFile[6*nu].replace(">","").rstrip()
        if dbFile[6*nu+1].rstrip() == "no site":
                return
        mutLociList = [i for i in dbFile[6*nu+1].rstrip().split()]
        nodePre = dbFile[6*nu+3].rstrip().split("-")
        mutationsList = [i for i in dbFile[6*nu+5].rstrip()]
        if loc in mutLociList:
                index = mutLociList.index(loc)
                if mutationsList[index] == mutation:
                        if rank == 1:
                                mut1EmergeNode[nodeName] = nodePre
                        else:
                                mut2EmergeNode[nodeName] = nodePre
                else:
                        pass
        else:
                pass

def mutationsOrder(mutation1Loc,mutation2Loc):
        mut1Node = mut1EmergeNode.keys()
        mut2Node = mut2EmergeNode.keys()
        mut12 = judgeOrder(mut1Node,mut2Node,rank=1)
        mut21 = judgeOrder(mut2Node,mut1Node,rank=2)    
        return mut12,mut21

def judgeOrder(mut1Node,mut2Node,rank=1):
        times = 0
        for node1 in mut1Node:
                for node2 in mut2Node:
                        if rank == 1:
                                if node1 in mut2EmergeNode.get(node2):
                                        times += 1
                                        if outNodeDict.get(node2):
                                                output1.write("%s\t%s\n" %(node1,outNodeDict.get(node2)))               
                                        else:
                                                output1.write("%s\t%s\n" %(node1,node2))
                                else:
                                        pass
                        else:
                                if node1 in mut1EmergeNode.get(node2):
                                        times += 1
                                        if outNodeDict.get(node2):
                                                output2.write("%s\t%s\n" %(node1,outNodeDict.get(node2)))
                                        else:
                                                output2.write("%s\t%s\n" %(node1,node2))
                                else:
                                        pass
        return times

def anoLoci(nodeList,rank=1):
        for node in nodeList:
                if outNodeDict.get(node):
                        if rank == 1:
                                mut1InTermianl.write(outNodeDict.get(node)+"\n")
                        else:
                                mut2InTermianl.write(outNodeDict.get(node)+"\n")        
                else:
                        if rank == 1:
                                mut1InNode.write(node+"\n")
                        else:
                                mut2InNode.write(node+"\n")

for nu in range(nodeNum):
        detect_mutation(mutation1Loc,mutation1,nu,rank=1)
        detect_mutation(mutation2Loc,mutation2,nu,rank=2)

mut12,mut21 = mutationsOrder(mutation1Loc,mutation2Loc)

anoLoci(mut1EmergeNode.keys(),rank=1)
anoLoci(mut2EmergeNode.keys(),rank=2)
print("%s before %s: %s" %(mut1,mut2,mut12))
print("%s before %s: %s" %(mut2,mut1,mut21))
