# -*- coding:Latin-1 -*-
"""
open the file created by ASE_eqtl_validation, add the counts for het and hom cross
"""

#############
### ASE 2 ###
#############

import os
import csv
import time
from operator import itemgetter


start_time = time.time()


#############
### PATHS ###
#############


OutName = "Treated"
eqtl_folder = ".\Data"
cnt_folder=".\GeneCnts"

file = ".\Data\\" + OutName + "_parsed_CrossInfo.txt"
out_file_path = ".\" + OutName + "_parsed_CrossInfo_CntData.txt"

print("input file", file)
print("output file", out_file_path)
input("press enter to continue")

#################
### Main loop ###
#################

### Open eQTL file and store data
print("read eqtl file")
data_ns=[]
first_line=0
with open(file,"r") as input_file:
    data1=input_file.readlines()
    for line in data1:
        tmp_line=line.split("\t")
        tmp_line[-1] = tmp_line[-1].strip() #header is now an array, without the end of line
        if first_line == 0:
            data_header=tmp_line
            first_line = 1
        else:
            data_ns.append(tmp_line)

# sort data
data=sorted(data_ns, key=itemgetter(3))

# Files with gene cnts
cnt_files=os.listdir(cnt_folder)

# indexes
cnt_gene_index=0
cnt_total_index=7
eqtl_gene_index=3

# loop through the cnt files
print("start looping through count files")
start_time=time.time()

for C_file in cnt_files:
    C_path=cnt_folder+"\\"+C_file
    #information on file
    sample_name="c" + C_file.split("_line")[0] #name of the sample
    cross="c" + C_file.split("_")[0] #name of the cross
    #name of values in data header to use to get index
    cross_status=cross+"_status"
    cross_l1=cross+"_l1"
    cross_l2=cross+"_l2"
    #name of values to write to header
    cross_alt=sample_name+"_alt"
    cross_ref=sample_name+"_ref"
    cross_tot=sample_name+"_total"
    #add information to header
    data_header.append(cross_alt)
    data_header.append(cross_ref)
    data_header.append(cross_tot)
    #index
    eqtl_status_index=data_header.index(cross_status) #this index give the column from data with the status of the cross
    eqtl_l1_index=data_header.index(cross_l1) #this is the line 1
    eqtl_l2_index=data_header.index(cross_l2) #this is for line 2
    #variable for file treatment
    tmp_cnts=[]
    first_line=0
    #open and store the file to tmp_cnts
    print("open", C_file)
    with open(C_path,"r") as C_input:
        C_data=C_input.readlines()
        for gene in C_data:
            tmp_line=gene.split("\t")
            tmp_line[-1] = tmp_line[-1].strip() #header is now an array, without the end of line
            if first_line == 0:
                tmp_header=tmp_line
                first_line = 1
            else:
                tmp_cnts.append(tmp_line)
    #start looping on data
    eqtl_iterator=0
    C_iterator=0
    while eqtl_iterator < len(data):
        #get information on eqtl in the cross
        eqtl_status=data[eqtl_iterator][eqtl_status_index] #het or other
        eqtl_gene=data[eqtl_iterator][eqtl_gene_index] #the eqtl gene
        #loop in cnt file
        while C_iterator < len(tmp_cnts):
            C_gene=tmp_cnts[C_iterator][cnt_gene_index]
            #check gene
            if eqtl_gene == C_gene:
                tot_cnt=tmp_cnts[C_iterator][cnt_total_index]
                #if eqtl is not heterozygote
                if eqtl_status == "Het":
                    eqtl_l1=data[eqtl_iterator][eqtl_l1_index] #the first line
                    eqtl_l2=data[eqtl_iterator][eqtl_l2_index] #the second line
                    #get information on line in cross
                    l1_type=str(eqtl_l1).split(":")[1]
                    l1_name=str(eqtl_l1).split(":")[0]
                    l2_type=str(eqtl_l2).split(":")[1]
                    l2_name=str(eqtl_l2).split(":")[0]
                    #with information, get the alt line and ref line, get index in tmp_cnts, store the value
                    if l1_type == "alt" and l2_type == "ref":
                        alt_header=tmp_header.index(l1_name)
                        alt_val=tmp_cnts[C_iterator][alt_header]
                        ref_header=tmp_header.index(l2_name)
                        ref_val=tmp_cnts[C_iterator][ref_header]
                    elif l1_type == "ref" and l2_type == "alt":
                        alt_header=tmp_header.index(l2_name)
                        alt_val=tmp_cnts[C_iterator][alt_header]
                        ref_header=tmp_header.index(l1_name)
                        ref_val=tmp_cnts[C_iterator][ref_header]
                    #add data, first alternate, then ref, then total
                    data[eqtl_iterator].append(alt_val)
                    data[eqtl_iterator].append(ref_val)
                    data[eqtl_iterator].append(tot_cnt)
                    break
                elif eqtl_status == "Hom_ref" or eqtl_status == "Hom_alt" :
                    #get line names
                    eqtl_l1=data[eqtl_iterator][eqtl_l1_index] #the first line
                    eqtl_l2=data[eqtl_iterator][eqtl_l2_index] #the second line
                    #get information on line in cross
                    l1_name=str(eqtl_l1).split(":")[0]
                    l2_name=str(eqtl_l2).split(":")[0]
                    #get counts
                    l1_header=tmp_header.index(l1_name)
                    l1_val=tmp_cnts[C_iterator][l1_header]
                    l2_header=tmp_header.index(l2_name)
                    l2_val=tmp_cnts[C_iterator][l2_header]
                    #add data
                    data[eqtl_iterator].append(l1_val)
                    data[eqtl_iterator].append(l2_val)
                    data[eqtl_iterator].append(tot_cnt)
                    break
                else:
                    #add two different variables for alt and ref count
                    data[eqtl_iterator].append("Hom_sample")
                    data[eqtl_iterator].append("Hom_sample")
                    data[eqtl_iterator].append(tot_cnt)
                    break
            if C_iterator == (len(tmp_cnts)-1):
                print("gene not found")
                print(eqtl_gene)
                print(C_gene)
                input("enter")
            C_iterator+=1
        eqtl_iterator+=1

print("\n\n\n######################")
print("END parsing data")
print("--- %s seconds ---" % (time.time() - start_time))

with open(out_file_path,"w") as out_file:
    out_file.write("\t".join(data_header) + "\n")
    for line in data:
        out_file.write("\t".join(line) + "\n")

print("\n\n\n######################")
print("END script")
print("--- %s seconds ---" % (time.time() - start_time))
