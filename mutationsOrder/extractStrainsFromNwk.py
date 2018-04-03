import re

"""
pre-process
treeFile = open(sys.argv[1]).read()
simplified_tree=re.compile('\d+\.\d+').sub('',treeFile).replace(';','').replace('-','').replace('_','')
"""
def substractNode2Terminal(processedTree,nodeName):
        output = open("Node"+str(nodeName),"a")
        startIndex = processedTree.index("Node%s:" %(nodeName))
        bracketRightNum = 0
        endIndex = 0
        for nu in range(startIndex-1,-1,-1):
                if processedTree[nu] == ")":
                        bracketRightNum += 1
                elif processedTree[nu] == "(":
                        bracketRightNum -= 1
                else:
                        pass
                if bracketRightNum == 0:
                        endIndex = nu
                        break
                else:
                        pass
        extractPart = processedTree[endIndex:startIndex]
        strainList = re.compile('(\)Node\d+|\(|\)|:)').sub("",extractPart).split(",")
        for strain in strainList:
                output.write(strain+"\n")


def substractNode2Node(processedTree,startNode,endNode):
        output = open("Node"+str(startNode)+"toNode"+str(endNode),"a")
        startIndex = min(processedTree.index("Node"+str(startNode)+":"),processedTree.index("Node"+str(endNode)+":"))
        bracketRightNum = 0
        endIndex = 0
        for nu in range(startIndex-1,-1,-1):
                if processedTree[nu] == ")":
                        bracketRightNum += 1
                elif processedTree[nu] == "(":
                        bracketRightNum -= 1
                else:
                        pass
                if bracketRightNum == 0:
                        endIndex = nu
                        break
                else:
                        pass
        extractPart = processedTree[endIndex:startIndex]
        strainList = re.compile('(\)Node\d+|\(|\)|:)').sub("",extractPart).split(",")
        for strain in strainList:
                output.write(strain+"\n")
