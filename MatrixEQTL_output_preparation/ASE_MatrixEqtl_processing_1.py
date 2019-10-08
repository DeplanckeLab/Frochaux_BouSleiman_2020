# -*- coding:Latin-1 -*-
# prepare .cis file into .txt for further processing

###########
### ASE ###
###########

import os
import csv
import argparse
import shutil #for removing folder
import re #regex for changing DGRP with line
import time
import tkinter as tk
from tkinter import filedialog
import ntpath
import itertools
import operator #for sorting

##############
# START TIME #
##############
start_time = time.time()


#################
### Functions ###
#################

"""
def get_SNP(eqtl_entry):
    "return information about snp"
    work=eqtl_entry.split("\t")
    chrom=work[0].split(":")[0]
    SNP_all=work[0].split(":")[1]
    pos=SNP_all.split("_")[0]
    var=SNP_all.split("_")[1]
    ass_gene=work[1]
    return(chrom, pos, var, ass_gene)
"""

##################
### USER INPUT ###
##################

### Input file winow ###
root = tk.Tk()
root.withdraw()

file = filedialog.askopenfilename() # get the complete path including the file name
file_name_complete = ntpath.basename(file) # get the name of the file including expansion

# get file path
find_separator = list(file) # get the separator used in file
for elt in find_separator:
    if elt == "/":
        path_separator = elt
        break
    if elt == "\\":
        path_separator = elt
        break

file_extension=file_name_complete.split(".")[1] # get the type of the file
file_name=file_name_complete.split(".")[0]

find_path=file.split(path_separator) # split the name of the path according to the separator found
find_path.remove(file_name_complete)
path=path_separator.join(find_path)

print(file)
OutName = input("name of the ouptut file:")
#print(OutName)

############
### Main ###
############

### create output and error file

# errorFileName=path + path_separator + OutName + "_error.txt"
outputFileName=path + path_separator + OutName + "_parsed.txt"
#error = open(errorFileName, "w")
out_file = open(outputFileName, "w")

# store parsed data in new array
ParsedData=[]
eqtl_cnt = 0

with open(file, "r") as cis_file: #, open(errorFileName, "w") as err_file:
    eqtl=cis_file.readlines() #read eqtls one by one
    iterator_eqtl=0
    print("file loaded")
    print("--- %s seconds ---" % (time.time() - start_time))
    while iterator_eqtl < len(eqtl): #loop through eqtl
        if iterator_eqtl == 0: #header processing
            eqtl_header=eqtl[iterator_eqtl]
            tmp_header=eqtl_header.split("\t")
            tmp_header[-1] = tmp_header[-1].strip() #header is now an array, without the end of line
            tmp_header.append("chrom")
            tmp_header.append("position")
            tmp_header.append("variant")
            out_file.write(tmp_header[-3]+ "\t" + tmp_header[-2] + "\t" +  tmp_header[-1] + "\t" + "\t".join(tmp_header[1:len(tmp_header)-3]) + "\n")
            iterator_eqtl+=1
            print("header processed")
        else:
            eqtl_cnt+=1
            #process line like the header
            eqtl_line=eqtl[iterator_eqtl]
            tmp_eqtl=eqtl_line.split("\t")
            tmp_eqtl[-1]=tmp_eqtl[-1].strip()
            e_chrom=tmp_eqtl[0].split(":")[0]
            e_pos=tmp_eqtl[0].split(":")[1].split("_")[0]
            e_var=tmp_eqtl[0].split(":")[1].split("_")[1]
            tmp_eqtl.append(e_chrom)
            tmp_eqtl.append(int(e_pos))
            tmp_eqtl.append(e_var)
            ParsedData.append(tmp_eqtl)
            iterator_eqtl+=1
            if iterator_eqtl % 500 == 0:
                print(iterator_eqtl, "lines processed out of", len(eqtl))
                print("--- %s seconds ---" % (time.time() - start_time))

# order data
ParsedData.sort(key = operator.itemgetter(6, 7))

# write to output file
for line in ParsedData:
    out_file.write(line[-3]+ "\t" + str(line[-2]) + "\t" +  line[-1] + "\t" + "\t".join(line[1:len(line)-3]) + "\n")

out_file.close()

print("\n\n\n######################")
print("--- %s seconds ---" % (time.time() - start_time))


########################
### Arguments parser ###
########################

"""
parser=argparse.ArgumentParser()
parser.add_argument("--eqtl", required=True, help="gene")
parser.add_argument("--out", required=False, default="out")
parser.add_argument("--Directory", required=False, default="out")
parser.add_argument("--vcf", required=False, default="out")
args = parser.parse_args()



print("args")
print("args.eqtl")
"""
