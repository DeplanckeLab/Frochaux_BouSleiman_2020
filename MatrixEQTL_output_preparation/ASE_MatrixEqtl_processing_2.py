# -*- coding:Latin-1 -*-

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

#########################
### Files and Folders ###
#########################
"""
needed:
eqtl file as outputed by matrix eqtl
vcf file that contains all the variants
cross file: to determin which cross to use for each eqtl
"""

# 1. vcf folders and files
# need to split vcf first for memory issues, vcf file can be found at DGRP freez 2: http://dgrp2.gnets.ncsu.edu/data.html
vcf_file= #to add
vcf_folder = #to add
tmp_folder = vcf_folder+"\\tmp"

# 2. files with informations on cross
cross_file = ".\Data\Cross_references.txt"

# 3. eqtl files
OutName = "Interaction"
eqtl_folder = ".\Data"
eqtl_file = ".\Data\\" + OutName + "_parsed.txt"

script_output_file = eqtl_folder + "\\" + OutName + "_parsed_CrossInfo.txt"
error_output_file = eqtl_folder + "\\" + OutName + "_error.txt"

#################
### Functions ###
#################

def get_SNP(eqtl_entry):
    "return information about snp"
    work=eqtl_entry.split("\t")
    chrom=work[0].split(":")[0]
    SNP_all=work[0].split(":")[1]
    pos=SNP_all.split("_")[0]
    var=SNP_all.split("_")[1]
    ass_gene=work[1]
    return(chrom, pos, var, ass_gene)


    

#################
### Split vcf ###
#################
"""
split vcf into chromosome based vcf
"""

#unmark to run
"""
vcf_file="C:\\Users\\Michael\\Documents\\Deplancke\\data\\ntc BRB-seq\\vcf\\dgrp2.vcf"

if not os.path.exists(tmp_folder):
    os.makedirs(tmp_folder)
print("temporary folder", tmp_folder, "created")

print("start processing vcf")

vcf_header=[]
cur_chr="lol"
line_cnt=0
with open(vcf_file, "r") as v:
    for line in v:
        work=str(line).split("\t")
        start=work[0]
        if start[0:2] == "##":
            vcf_header.append(line)
        elif start[0:7] == "#CHROM":
            vcf_ref=line
        else:
            #1. check if there is a new chromosome, if yes:
            #       close current open file
            #       swith current chromosome
            #       open write file
            if start != cur_chr:
                if cur_chr == "lol":
                    #switch chromosome
                    print("start chromosome: ", start)
                    cur_chr=start
                    #create new file path
                    cur_fpath=tmp_folder+"\\chrom_"+cur_chr+".vcf"
                    output_file=open(cur_fpath, "w")
                    #write
                    output_file.write(vcf_ref)
                    output_file.write(line)
                else:
                    #switch chromosome
                    print("start chromosome: ", start)
                    cur_chr=start
                    #close current output file
                    output_file.close()
                    #create new file path
                    cur_fpath=tmp_folder+"\\chrom_"+cur_chr+".vcf"
                    output_file=open(cur_fpath, "w")
                    #write
                    output_file.write(vcf_ref)
                    output_file.write(line)
            else:
                output_file.write(line)
        if line_cnt % 100000 == 0:
            print("processed: ", line_cnt, "lines")
        line_cnt+=1



#shutil.rmtree(tmp_folder, ignore_errors=True)

output_file.close()
input("press enter to continue") 
"""

############
### Main ###
############


#list of files from the vcf processing, used if vcf file is processed
vcf_list = os.listdir(tmp_folder)

### parsing cross file
i=0
cross_ref=[]
with open(cross_file, "r") as c_file:
    cross=c_file.readlines()
    for cr in cross:
        if i == 0:
            cross_header=cr
        else:
            cross_ref.append(cr.split("\t"))
        i=1
### end parsing cross file

i=0 #to identify first line
cur_chrom="lol" #initiating first file
eqtl_cnt=0
nb_found=0
nb_NotFound=0

with open(eqtl_file, "r") as eq_file, open(script_output_file,"w") as out_file, open(error_output_file, "w") as err_file:
    eqtl=eq_file.readlines() #read eqtls one by one
    iterator_eqtl=0
    print("eqtl file loaded")
    print(eqtl_file)
    print("--- %s seconds ---" % (time.time() - start_time))
    while iterator_eqtl < len(eqtl): #loop through eqtl
        if iterator_eqtl == 0: #header processing
            eqtl_header=eqtl[iterator_eqtl]
            tmp_header=eqtl_header.split("\t")
            tmp_header[-1] = tmp_header[-1].strip() #header is now an array, without the end of line
            err_file.write(tmp_header[0] + "\t" + "\t".join(tmp_header[1:]) + "\tstatus\n")
            for cr in cross_ref:
                tmp_header.append("c"+cr[0]+"_status") #to header, add info about cross status
                tmp_header.append("c"+cr[0]+"_l1") #line1
                tmp_header.append("c"+cr[0]+"_l2") #line2
            out_file.write(tmp_header[0] + "\t" + "\t".join(tmp_header[1:]) + "\n") #write to output file
            print("eatl header processed")
            iterator_eqtl+=1
        else:
            eqtl_cnt+=1
            #process line like the header
            eqtl_line=eqtl[iterator_eqtl]
            tmp_eqtl=eqtl_line.split("\t")
            tmp_eqtl[-1]=tmp_eqtl[-1].strip()
            e_chrom=tmp_eqtl[0]
            e_pos=tmp_eqtl[1]
            e_var=tmp_eqtl[2]
            e_assGene=tmp_eqtl[3]
            if iterator_eqtl % 500 == 0:
                print(iterator_eqtl, "eqtl processed out of", len(eqtl))
                print("--- %s seconds ---" % (time.time() - start_time))
            #if new chromosome, open new file and store into current vcf
            if e_chrom != cur_chrom:
                iterator_vcf=0
                print("change chromosome, from", cur_chrom, "to", e_chrom)
                cur_chrom=e_chrom
                vcf_file = open(tmp_folder+"\\chrom_"+cur_chrom+".vcf", "r")
                vcf = vcf_file.readlines()
                vcf_file.close
                """
                with open(tmp_folder+"\\chrom_"+cur_chrom+".vcf", "r") as vcf_file:
                    vcf=vcf_file.readlines()
                """
                print(cur_chrom, "loaded")
            #browse through variant in the vcf file
            found=0
            while iterator_vcf < len(vcf):
                if found:
                    iterator_eqtl+=1
                    break
                #remove header
                if iterator_vcf == 0:
                    vcf_header=vcf[iterator_vcf]
                    tmp_vcf=vcf_header.split("\t")
                    tmp_vcf[-1] = tmp_vcf[-1].strip()
                    print("header processed")
                else:
                    work=vcf[iterator_vcf].split("\t")
                    if int(e_pos) < int(work[1]):
                        err_file.write(tmp_eqtl[0] + "\t" + "\t".join(tmp_eqtl[1:]) + "\tnot found\n")
                        tmp_out="NA"
                        for cr in cross_ref:
                            tmp_eqtl.append(tmp_out) #fill the 3 fields
                            tmp_eqtl.append(tmp_out)
                            tmp_eqtl.append(tmp_out)
                        out_file.write(tmp_eqtl[0] + "\t" + "\t".join(tmp_eqtl[1:]) + "\n")
                        found=True
                        nb_NotFound+=1
                    #check for position
                    if int(work[1]) == int(e_pos):
                        v_var=work[3]+work[4]
                        #check if same variant
                        if v_var == e_var.replace("/", ""):
                            err_file.write(tmp_eqtl[0]+ "\t" + "\t".join(tmp_eqtl[1:]) + "\tfound\n")
                            for cr in cross_ref:
                                #for each cross, get the lines and the variants
                                #1/1 is alternate, 0/0 is reference
                                l1 = cr[1]
                                l1_L=cr[1].replace("DGRP-", "line_")
                                l2 = cr[2]
                                l2_L=cr[2].replace("DGRP-", "line_")
                                #get the index of lines
                                l1_index=tmp_vcf.index(l1_L)
                                l2_index=tmp_vcf.index(l2_L)
                                #get the variants for each line
                                l1_var = work[l1_index]
                                l2_var = work[l2_index]
                                #change output type
                                l1_out="NA"
                                l2_out="NA"
                                if l1_var == "1/1":
                                    l1_out = "alt"
                                if l1_var == "0/0":
                                    l1_out = "ref"
                                if l2_var == "1/1":
                                    l2_out = "alt"
                                if l2_var == "0/0":
                                    l2_out = "ref"
                                #change cross type
                                c_type=""
                                if l1_var == l2_var:
                                    c_type="Hom"+"_"+l1_out
                                if l1_var != l2_var:
                                    c_type="Het"
                                if l1_var == "./." or l2_var == "./.":
                                    c_type="NU"
                                tmp_eqtl.append(c_type)
                                tmp_eqtl.append(l1_L+":"+l1_out)
                                tmp_eqtl.append(l2_L+":"+l2_out)  
                            out_file.write(tmp_eqtl[0] + "\t" + "\t".join(tmp_eqtl[1:]) + "\n")
                            nb_found+=1
                            iterator_vcf-=1
                            found=True
                iterator_vcf+=1


print("\n\n\n######################")
print("Stats:\nTotal:", len(eqtl),"\nFound:", nb_found, "\nNot Found", nb_NotFound)
print("--- %s seconds ---" % (time.time() - start_time))
