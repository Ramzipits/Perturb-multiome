#!/bin/python3

from os import listdir
import os
import subprocess
import gzip

file_name_list = [f_name for f_name in listdir("/lustre2/sctxt_k2") if ".txt" in f_name]
wdir = "/lustre2/sctxt_k2"

guide_list = [guide.strip("\n") for guide in open("Z_37_GBC_count.txt")]
print guide_list


counts_dic = {}
for i in range(0,37):
	counts_dic[i] = []

for file_name in file_name_list:
	with open(file_name,'r') as f:
		print file_name
		for guide_id in range(0,37):
			cmd_fwd = "cat " + file_name + "| grep " + guide_list[guide_id] + " | sort | uniq | wc -l "
			ps = subprocess.Popen(cmd_fwd, shell=True, stdout=subprocess.PIPE)
			output_fwd = int(ps.communicate()[0].strip("\n"))
			counts_dic[guide_id].append(output_fwd)

outf = open("GBClib-allcounts.txt","w")

for file_id in range(len(file_name_list)):
	print file_id 
	name = file_name_list[file_id]
	outf.write(name.split(".")[0] + "\t") 
	for key in counts_dic.keys():
		outf.write(str(counts_dic[key][file_id]) + "\t")
	outf.write("\n")


outf.close()



