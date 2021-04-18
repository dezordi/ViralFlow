#!/usr/bin/python3
# -*- coding: utf-8 -*-

__author__ = "Filipe Z. Dezordi"
__license__ = "GPL"
__version__ = "0.0.1"
__maintainer__ = "Filipe Z. Dezordi"
__email__ = "zimmer.filipe@gmail.com"
__date__ = "2021/11/11"
__username__ = "dezordi"

import argparse, subprocess, shlex, csv, re, os
import pandas as pd
import numpy as np

parser = argparse.ArgumentParser(description = 'This script returns de major and minor variants per positions')
parser.add_argument("-in","--input_file",help="TSV minor formated file", required=True)
args = parser.parse_args()
input_file = args.input_file
output_tsv = open(input_file+'.fmt','w')
output_tsv_writer = csv.writer(output_tsv, delimiter='\t')
output_tsv_writer.writerow(["GENOME","POS","REGION","DEPTH","A_DEPTH","A_PLUS","A_MINUS","C_DEPTH","C_PLUS","C_MINUS","G_DEPTH","G_PLUS","G_MINUS","T_DEPTH","T_PLUS","T_MINUS","INDEL_DEPTH","INDEL_PLUS","INDEL_MINUS","INDEL","PUTATIVE_MINOR","MAJOR","MINOR","MAJOR_DEPTH","MINOR_DEPTH"])
out_list_lists = []
with open(input_file,'r') as minor_tsv:
    major = ''
    minor = ''
    major_depth = ''
    minor_depth = ''
    minor_tsv_reader = csv.reader(minor_tsv, delimiter='\t')
    for line in minor_tsv_reader:
        if line[20] == "A/C": 
            if int(line[4]) > int(line[7]):
                major = 'A'
                minor = 'C'
                major_depth = line[4]
                minor_depth = line[7]
            else:
                major = 'C'
                minor = 'A'
                major_depth = line[7]
                minor_depth = line[4]
        elif line[20] == "A/G":
            if int(line[4]) > int(line[10]):
                major = 'A'
                minor = 'G'
                major_depth = line[4]
                minor_depth = line[10]
            else:
                major = 'G'
                minor = 'A'
                major_depth = line[10]
                minor_depth = line[4]
        elif line[20] == "A/T":
            if int(line[4]) > int(line[13]):
                major = 'A'
                minor = 'T'
                major_depth = line[4]
                minor_depth = line[13]
            else:
                major = 'T'
                minor = 'A'
                major_depth = line[13]
                minor_depth = line[4]
        elif line[20] == "A/INDEL":
            if '-' in line[19]:
                if int(line[4]) > int(line[16]):
                    major = 'A'
                    minor = '-'
                    major_depth = line[4]
                    minor_depth = line[16]
                else:
                    major = '-'
                    minor = 'A'
                    major_depth = line[16]
                    minor_depth = line[4]
            else:
                if int(line[4]) > int(line[16]):
                    major = 'A'
                    minor = 'A'+line[19][1:]
                    major_depth = line[4]
                    minor_depth = line[16]
                else:
                    major = 'A'+line[19][1:]
                    minor = 'A'
                    major_depth = line[16]
                    minor_depth = line[4]
        elif line[20] == "T/INDEL":
            if '-' in line[19]:
                if int(line[13]) > int(line[16]):
                    major = 'T'
                    minor = '-'
                    major_depth = line[13]
                    minor_depth = line[16]
                else:
                    major = '-'
                    minor = 'T'
                    major_depth = line[16]
                    minor_depth = line[13]
            else:
                if int(line[13]) > int(line[16]):
                    major = 'T'
                    minor = 'T'+line[19][1:]
                    major_depth = line[13]
                    minor_depth = line[16]
                else:
                    major = 'T'+line[19][1:]
                    minor = 'T'
                    major_depth = line[16]
                    minor_depth = line[13]
        elif line[20] == "C/G":
            if int(line[7]) > int(line[10]):
                major = 'C'
                minor = 'G'
                major_depth = line[7]
                minor_depth = line[10]
            else:
                major = 'G'
                minor = 'C'
                major_depth = line[10]
                minor_depth = line[7]
        elif line[20] == "T/C":
            if int(line[13]) > int(line[7]):
                major = 'T'
                minor = 'C'
                major_depth = line[13]
                minor_depth = line[7]
            else:
                major = 'C'
                minor = 'T'
                major_depth = line[7]
                minor_depth = line[13]
        elif line[20] == "T/G":
            if int(line[13]) > int(line[10]):
                major = 'T'
                minor = 'G'
                major_depth = line[13]
                minor_depth = line[10]
            else:
                major = 'G'
                minor = 'T'
                major_depth = line[10]
                minor_depth = line[13]
        elif line[20] == "C/INDEL":
            if '-' in line[19]:
                if int(line[7]) > int(line[16]):
                    major = 'C'
                    minor = '-'
                    major_depth = line[7]
                    minor_depth = line[16]
                else:
                    major = '-'
                    minor = 'C'
                    major_depth = line[16]
                    minor_depth = line[7]
            else:
                if int(line[7]) > int(line[16]):
                    major = 'C'
                    minor = 'C'+line[19][1:]
                    major_depth = line[7]
                    minor_depth = line[16]
                else:
                    major = 'C'+line[19][1:]
                    minor = 'C'
                    major_depth = line[16]
                    minor_depth = line[7]

        elif line[20] == "G/INDEL":
            if '-' in line[19]:
                if int(line[10]) > int(line[16]):
                    major = 'G'
                    minor = '-'
                    major_depth = line[10]
                    minor_depth = line[16]
                else:
                    major = '-'
                    minor = 'G'
                    major_depth = line[16]
                    minor_depth = line[10]
            else:
                if int(line[10]) > int(line[16]):
                    major = 'G'
                    minor = 'G'+line[19][1:]
                    major_depth = line[10]
                    minor_depth = line[16]
                else:
                    major = 'G'+line[19][1:]
                    minor = 'G'
                    major_depth = line[16]
                    minor_depth = line[10]
        else:
            major = 'CURATION'
            minor = 'CURATION'
            major_depth = ''
            minor_depth = ''
        if 'GENOME' in line[0]:
            pass
        else:
            out_list_lists.append([line[0],line[1],line[2],line[3],line[4],line[5],line[6],line[7],line[8],line[9],line[10],line[11],line[12],line[13],line[14],line[15],line[16],line[17],line[18],line[19],line[20],major,minor,major_depth,minor_depth])
output_tsv_writer.writerows(out_list_lists)