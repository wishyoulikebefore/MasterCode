import sys
from collections import defaultdict

strainList=open(sys.argv[1]).readlines()
output=open("resistantStrains_mutations_info","a")
DR_list=open("DR.list","a")
DR_RIF_list=open("DR_RIF.list","a")
DR_INH_list=open("DR_INH.list","a")
DR_FLUO_list=open("DR_FLUO.list","a")
DR_PZA_list=open("DR_PZA.list","a")
DR_EMB_list=open("DR_EMB.list","a")
MDR_list=open("MDR.list","a")
MDR_PZA_list=open("MDR_PZA.list","a")
MDR_EMB_list=open("MDR_EMB.list","a")
pre_XDR_list=open("pre-XDR.list","a")
pre_XDR_PZA_list=open("pre-XDR_PZA.list","a")
pre_XDR_EMB_list=open("pre-XDR_EMB.list","a")
XDR_list=open("XDR.list","a")
XDR_PZA_list=open("XDR_PZA.list","a")
XDR_EMB_list=open("XDR_EMB.list","a")
S_list=open("S.list","a")
INH_mutation_spectrum = open("INH_mutation_spectrum","a")
RIF_mutation_spectrum = open("RIF_mutation_spectrum","a")
EMB_mutation_spectrum = open("EMB_mutation_spectrum","a")
PZA_mutation_spectrum = open("PZA_mutation_spectrum","a")

def MutationProcess(strain):
    resistantFile = open(strain+".txt").readlines()
    mutationDict = defaultdict(set)
    for row in resistantFile:
        rowOfList = row.rstrip().split()
        mutation = rowOfList[1]+"_"+rowOfList[5]
        if "RIF" in row:
            mutationDict["RIF"].add(mutation)
        elif "ISO" in row:
            mutationDict["INH"].add(mutation)
        elif "FLUOROQUINOLONES" in row:
            mutationDict["FLUO"].add(mutation)
        elif "AMIKACIN" in row or "CAPREOMYCIN" in row or "KANAMYCIN" in row:
            mutationDict["ACK"].add(mutation)
        elif "ETHAMBUTOL" in row:
            mutationDict["EMB"].add(mutation)
        elif "PYRAZINAMIDE" in row:
            mutationDict["PZA"].add(mutation)
        else:
            pass
    RIF = ["+".join(mutationDict.get("RIF")) if mutationDict.get("RIF") else ""]
    INH = ["+".join(mutationDict.get("INH")) if mutationDict.get("INH") else ""]
    ACK = ["+".join(mutationDict.get("ACK")) if mutationDict.get("ACK") else ""]
    FLUO = ["+".join(mutationDict.get("FLUO")) if mutationDict.get("FLUO") else ""]
    EMB = ["+".join(mutationDict.get("EMB")) if mutationDict.get("EMB") else ""]
    PZA = ["+".join(mutationDict.get("PZA")) if mutationDict.get("PZA") else ""]
    period = judgeResistancePeriod(INH,RIF,ACK,FLUO,EMB,PZA)
    output.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %(strain,INH,RIF,ACK,FLUO,EMB,PZA,period))
    distribute(strain,period,INH,RIF,ACK,FLUO,EMB,PZA)
    IREP_calculate(strain,period,INH,RIF,EMB,PZA)
    print("%s.txt has been processed" % (strain))

def judgeResistancePeriod(INH,RIF,ACK,FLUO,EMB,PZA):
    if RIF != "" and INH != "":
        if FLUO == "" and ACK == "":
            return "MDR"
        elif FLUO != "" and ACK != "":
            return "XDR"
        else:
            return "pre-XDR"
    elif FLUO == "" and ACK == "" and INH == "" and RIF == "" and PZA == "" and EMB == "":
        return "S"
    else:
        return "DR"

def distribute(strain,period,INH,RIF,ACK,FLUO,EMB,PZA):
    if period == "MDR":
        MDR_list.write(strain + "\n")
        if PZA != "":
            MDR_PZA_list.write(strain + "\n")
        if EMB != "":
            MDR_EMB_list.write(strain + "\n")
    elif period == "pre-XDR":
        pre_XDR_list.write(strain + "\n")
        if PZA != "":
            pre_XDR_PZA_list.write(strain + "\n")
        if EMB != "":
            pre_XDR_EMB_list.write(strain + "\n")
    elif period == "XDR":
        XDR_list.write(strain + "\n")
        if PZA != "":
            XDR_PZA_list.write(strain + "\n")
        if EMB != "":
            XDR_EMB_list.write(strain + "\n")
    elif period == "DR":
        DR_list.write(strain + "\n")
        if RIF != "":
            DR_RIF_list.write(strain + "\n")
        if INH != "":
            DR_INH_list.write(strain + "\n")
        if PZA != "":
            DR_PZA_list.write(strain + "\n")
        if EMB != "":
            DR_EMB_list.write(strain + "\n")
    else:
        S_list.write(strain + "\n")

def IREP_calculate(strain,period,INH,RIF,EMB,PZA):
    if INH != "":
        INH_mutation_spectrum.write("%s\t%s\t%s\n" %(strain,period,INH))
    if RIF != "":
        RIF_mutation_spectrum.write("%s\t%s\t%s\n" %(strain,period,RIF))
    if EMB != "":
        EMB_mutation_spectrum.write("%s\t%s\t%s\n" %(strain,period,EMB))
    if PZA != "":
        PZA_mutation_spectrum.write("%s\t%s\t%s\n" %(strain,period,PZA))

if __name__ == "main":
    for strain in strainList:
        strain = strain.rstrip()
        MutationProcess(strain)
