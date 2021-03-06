import sys
import re
from collections import Iterable

tree=""
treeFile = open(sys.argv[1]).readlines()
simplified_tree=re.compile('(\d+\.\d+|E\-\d)').sub('',tree).replace(';','').replace('-','').replace('_','').replace(':','')
strains_input=[]
for row in open(sys.argv[2]).readlines():
        strains_input.append(row.rstrip().replace('-','').replace('_',''))
strains_in_tree=simplified_tree.replace('(','').replace(')','')[19:].split(',')

def ei(strains):
        if not isinstance(strains,list):
                return "Input error"
        if len(strains)<2:
                return "Input strains < 2"
        name_len=0
        name_sum={}
        for row in strains:
                name_len+=len(row)
                name_sum[row]=simplified_tree.index(row)
        sort_name=sorted(name_sum.items(),key=lambda e:e[1])
        min_index=sort_name[0][1]-1
        max_index=sort_name[-1][1]+len(sort_name[-1][0])
        if simplified_tree[min_index] == "(" and simplified_tree[max_index] == ")" :
                substract=''
                for i in range(min_index,max_index+1):
                        substract+=simplified_tree[i] 
                if substract.count("(") != substract.count(")"):
                        if substract.count("(") > substract.count(")") and simplified_tree[max_index+1] != ")" :
                                return "Evovle independently !"
                        if substract.count("(") < substract.count(")") and simplified_tree[min_index-1] != "(" :
                                return "Evovle independently !"
                if len(substract.replace(',','').replace('(','').replace(')','')) == name_len:
                        return "Not evovle independently !"
                else:
                        return "Evovle independently !"
        else:
                return "Evovle independently !"

print(sys.argv[2]+ " " +ei(strains_input))

def positive(adjacent_strains):
        count=1
        i=2
        while i <= len(adjacent_strains):
                if ei(adjacent_strains[:2]) != "Evovle independently !":
                                        i+=1
                else:
                        count+=1
                        adjacent_strains=adjacent_strains[i-1:]
                        i=2
        return count

def negative(adjacent_strains):
        count=1
        j=len(adjacent_strains)
        while j >=2 :
                if ei(adjacent_strains[j-2:]) != "Evovle independently !":
                        j-=1
                else:
                        count+=1
                        adjacent_strains=adjacent_strains[:j-1]
                        j=len(adjacent_strains)
        return count

def adjacent_judge(adjacent_strains):
        if ei(adjacent_strains) != "Evovle independently !":
                return 1
        else:
                return min(negative(adjacent_strains),positive(adjacent_strains))                       
        
if ei(strains_input) == "Evovle independently !" :
        count=0
        sort_index=sorted([strains_in_tree.index(i) for i in strains_input])
        adjacent_substract=[sort_index[i]-sort_index[i-1] for i in range(1,len(sort_index))]+[0]
        separate=[]
        adjacent_strains=[]
        for i in range(len(adjacent_substract)):
                if adjacent_substract[i]==1:
                        adjacent_strains.append(strains_in_tree[sort_index[i]])
                elif adjacent_substract[i-1]==1:
                        adjacent_strains.append(strains_in_tree[sort_index[i]])
                        separate.append(adjacent_strains)
                        adjacent_strains=[]
                else:
                        adjacent_strains.append(strains_in_tree[sort_index[i]])
                        separate.append(adjacent_strains)
                        adjacent_strains=[]
        for item in separate:
                if len(item) == 1:
                        count+=1
                else:
                        count+=adjacent_judge(item)             
        print(sys.argv[2]+ " Evovle independently " + str(count) + " times !")
