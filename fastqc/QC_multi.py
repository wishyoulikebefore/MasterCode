#usage  p3 .py list -R/-noR
import sys
from collections import OrderedDict

output=open("QC_summary","a")
output.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %("file","Total_reads","Sequence_length","%GC","baseLtQ20_location","baseLtQ20_proportion","sequence_Q20","sequence_Q30","Fail"))

f1=open(sys.argv[1]).readlines()
f2=sys.argv[2]
strain_list=[i.rstrip() for i in f1]
summary_file="summary.txt"
data_file="fastqc_data.txt"

def count_sequence_quality(quality_dict):
	base_num = 0
	Q20 = 0
	Q30 = 0
	for key in quality_dict:
		if int(key) < 20:
			Q20 += int(quality_dict[key])
		if int(key) < 30:
			Q30 += int(quality_dict[key])
		base_num += int(quality_dict[key])
	Q20 = 1-Q20/base_num
	Q30 = 1-Q30/base_num
	return "%.4f"%(Q20),"%.4f"%(Q30)

def count_base_quality(quality_dict):
	location_list = []
	start = 0
	end = 0
	Q20_length = 0
	for key in quality_dict:
		if "-" not in key:
			location_list.append(key)
			Q20_length += 1
		else:
			keyOflist = key.split("-")
			Q20_length += int(keyOflist[1]) - int(keyOflist[0])
			if start == 0:
				start = int(keyOflist[0])
				end = int(keyOflist[1])
			else:
				if int(keyOflist[0]) == end+1:
					end = int(keyOflist[1])
				else:
					location_list.append("%s-%s"%(start,end))
					start = int(keyOflist[0])
					end = int(keyOflist[1])
	if len(location_list) == 0:
		return ("%s-%s"%(start,end)), Q20_length
	else:
		location_list.append("%s-%s"%(start,end))
		return location_list, Q20_length
	
def substract_data(file_path):
	F = open(file_path).readlines()
	row_num = len(F)
	page_end = 0
	base_quality_dict = OrderedDict()
	sequence_quality_dict = {}
	for nu in range(row_num+1):
		if F[nu].startswith("Total Sequences"):
			Total_reads = F[nu].rstrip().split()[2]
			page_end = nu
			break
		else:
			pass
	for nu in range(page_end+1,row_num+1):
		if F[nu].startswith("Sequence length"):
			Sequence_length = F[nu].rstrip().split()[2]
			if "-" in F[nu]:
				max_length = int(F[nu].rstrip().split("-")[1])
				page_end = nu
				break
			else:
				max_length = int(F[nu].rstrip().split()[-1])
				page_end = nu
				break
		else:
			pass
	for nu in range(page_end+1,row_num+1):
		if F[nu].startswith("%GC"):
			GC_content = F[nu].rstrip().split()[1]
			page_end = nu
			break
		else:
			pass
	for nu in range(page_end+1,row_num+1):
		if F[nu].startswith("#Base"):
			for num in range(nu+1,row_num+1):
				if F[num].startswith(">>END_MODULE"):
					page_end = num-1
					break
				else:
					rowOflist = F[num].rstrip().split()
					if round(float(rowOflist[1]),0) < 20:
						base_quality_dict[rowOflist[0]] = rowOflist[1]
					else:
						pass
			break
		else:
			pass
	if len(base_quality_dict) == 0:
		base_Q20_location = "All>20"
		base_Q20_length = 0
	else:
		base_Q20_location, base_Q20_length = count_base_quality(base_quality_dict)  
	for nu in range(page_end+1,row_num+1):
		if F[nu].startswith("#Quality"):
			for num in range(nu+1,row_num+1):
				if F[num].startswith(">>END_MODULE"):
					break
				else:
					row = F[num].rstrip().split()
					countNum = row[1].replace(".0","")
					if "E" in countNum:
						baseNum = countNum[:-2]
						fold = countNum[-1]
						countNum = float(baseNum)*(10**int(fold))
					sequence_quality_dict[row[0]] = countNum
			break
		else:
			pass
	Q20,Q30 = count_sequence_quality(sequence_quality_dict)
	base_Q20_proportion = "%.2f" %(base_Q20_length/max_length)
	return Total_reads,Sequence_length,GC_content,base_Q20_location,base_Q20_proportion,Q20,Q30

def FAIL_PASS(file_path):
	fail_index = []
	for row in open(file_path).readlines():
		row = row.rstrip().split("\t")
		if row[0] == "FAIL":
			fail_index.append(row[1])
		else:
			pass
	return str(fail_index)

def write_out(file_path):
	file1 = file_path+"/"+data_file
	file2 = file_path+"/"+summary_file
	output.write(file_path+"\t")
	output.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t" %(substract_data(file1)))
	output.write("%s\n" %(str(FAIL_PASS(file2))))
					
for strain in strain_list:
	if f2 == "R":
		path1 = strain+"_R1_fastqc"
		write_out(path1)
		path2 = strain+"_R2_fastqc"
		write_out(path2)
	else:
		path1 = strain+"_1_fastqc"
		write_out(path1)
		path2 = strain+"_2_fastqc"
		write_out(path1)

