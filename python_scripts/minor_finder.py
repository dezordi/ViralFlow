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

parser = argparse.ArgumentParser(description = 'This scripts returns putative minor variants, based on nucleotide diversity of reads mapped against a reference genome')
parser.add_argument("-in","--bam_readcount",help="bam-readcount output file", required=True)
parser.add_argument("-lm","--per_limit",help="Percentage limit of Depth Base / Total Base to consider as putative minor variant, default = 0.05, which represents 5 percent", type=float,default= 0.05)
parser.add_argument("-dp","--depth_number",help="Minimum depth value to consider as minor variant, default = 100", type=int, default= 100)
args = parser.parse_args()
#store arguments
bam_rc_file = args.bam_readcount
per_limit = args.per_limit
depth_value = args.depth_number
#create prefix
prefix = bam_rc_file
prefix = re.sub(".*\/","",prefix)
prefix = re.sub("\..*","",prefix)
with open(bam_rc_file,'r') as bc_file, open(bam_rc_file+'.fmt.tsv','w') as bc_formated_output:
    output_csv_writer = csv.writer(bc_formated_output,delimiter='\t')
    output_csv_writer.writerow(['GENOME','POS','REGION','DEPTH','A_DEPTH','A_PLUS','A_MINUS','C_DEPTH','C_PLUS','C_MINUS','G_DEPTH','G_PLUS','G_MINUS','T_DEPTH','T_PLUS','T_MINUS','INDEL_DEPTH','INDEL_PLUS','INDEL_MINUS','INDEL'])
    out_list = []
    del_dict = {'pos':[],'deletion':[],'del_depth':[],'del_depth_plus':[],'del_depth_minus':[]}
    bam_rc_file_read = csv.reader(bc_file, delimiter='\t')
    for line in bam_rc_file_read:
        var_pos = line[1].rstrip('\n')
        var_depth = line[3].rstrip('\n')
        A_depth_line = line[5].rstrip('\n')
        A_depth = ":".join(A_depth_line.split(":", 2)[:2])
        A_depth = re.sub(r".*:","",A_depth)
        A_plus = ":".join(A_depth_line.split(":", 6)[:6])
        A_plus = re.sub(r".*:","",A_plus)
        A_minus = ":".join(A_depth_line.split(":", 7)[:7])
        A_minus = re.sub(r".*:","",A_minus)
        C_depth_line = line[6].rstrip('\n')
        C_depth = ":".join(C_depth_line.split(":", 2)[:2])
        C_depth = re.sub(r".*:","",C_depth)
        C_plus = ":".join(C_depth_line.split(":", 6)[:6])
        C_plus = re.sub(r".*:","",C_plus)
        C_minus = ":".join(C_depth_line.split(":", 7)[:7])
        C_minus = re.sub(r".*:","",C_minus)
        G_depth_line = line[7].rstrip('\n')
        G_depth = ":".join(G_depth_line.split(":", 2)[:2])
        G_depth = re.sub(r".*:","",G_depth)
        G_plus = ":".join(G_depth_line.split(":", 6)[:6])
        G_plus = re.sub(r".*:","",G_plus)
        G_minus = ":".join(G_depth_line.split(":", 7)[:7])
        G_minus = re.sub(r".*:","",G_minus)
        T_depth_line = line[8].rstrip('\n')
        T_depth = ":".join(T_depth_line.split(":", 2)[:2])
        T_depth = re.sub(r".*:","",T_depth)
        T_plus = ":".join(T_depth_line.split(":", 6)[:6])
        T_plus = re.sub(r".*:","",T_plus)
        T_minus = ":".join(T_depth_line.split(":", 7)[:7])
        T_minus = re.sub(r".*:","",T_minus)
        major_depth_column = 0
        major_depth = 0
        if len(line) > 11: #bamreadcount output can have 11 or more columns, depending of number of putative indels
            for i in range(10,len(line)):
                #print(line[1],line[i])
                indel_depth_line = line[i].rstrip('\n')
                indel_seq = re.sub(r":.*","",indel_depth_line)
                indel_depth = ":".join(indel_depth_line.split(":", 2)[:2])
                indel_depth = re.sub(r".*:","",indel_depth)
                if int(indel_depth) > int(major_depth):
                    major_depth = indel_depth
                    major_depth_column = i
            indel_depth_line = line[major_depth_column].rstrip('\n')
            indel_seq = re.sub(r":.*","",indel_depth_line)
            indel_depth = ":".join(indel_depth_line.split(":", 2)[:2])
            indel_depth = re.sub(r".*:","",indel_depth)
            indel_plus = ":".join(indel_depth_line.split(":", 6)[:6])
            indel_plus = re.sub(r".*:","",indel_plus)
            indel_minus = ":".join(indel_depth_line.split(":", 7)[:7])
            indel_minus = re.sub(r".*:","",indel_minus)
            if '-' in indel_seq:
                del_length = len(indel_seq)-1
                #print(indel_seq,del_length,line[1])
                for i in range(0,del_length):
                    del_dict['pos'].append(int(line[1])+i)
                    del_dict['deletion'].append(indel_seq)
                    del_dict['del_depth'].append(indel_depth)
                    del_dict['del_depth_plus'].append(indel_plus)
                    del_dict['del_depth_minus'].append(indel_minus)
            #print(line[1],indel_depth,indel_plus,indel_minus,indel_seq)
        elif len(line) == 11:
            indel_depth_line = line[10].rstrip('\n')
            indel_seq = re.sub(r":.*","",indel_depth_line)
            indel_depth = ":".join(indel_depth_line.split(":", 2)[:2])
            indel_depth = re.sub(r".*:","",indel_depth)
            indel_plus = ":".join(indel_depth_line.split(":", 6)[:6])
            indel_plus = re.sub(r".*:","",indel_plus)
            indel_minus = ":".join(indel_depth_line.split(":", 7)[:7])
            indel_minus = re.sub(r".*:","",indel_minus)
            if '-' in indel_seq:
                del_length = len(indel_seq)-1
                #print(indel_seq,del_length,line[1])
                for i in range(0,del_length):
                    del_dict['pos'].append(int(line[1])+i)
                    del_dict['deletion'].append(indel_seq)
                    del_dict['del_depth'].append(indel_depth)
                    del_dict['del_depth_plus'].append(indel_plus)
                    del_dict['del_depth_minus'].append(indel_minus)
            #print(line[1],indel_depth,indel_plus,indel_minus,indel_seq)
        if int(var_pos) <= 265:
            region = '5UTR'
        elif int(var_pos) <= 21555:
            region = 'ORF1AB'
        elif int(var_pos) <= 25384:
            region = 'S'
        elif int(var_pos) <= 26220:
            region = 'ORF3A'
        elif int(var_pos) <= 26472:
            region = 'E'
        elif int(var_pos) <= 27191:
            region = 'M'
        elif int(var_pos) <= 27387:
            region = 'ORF6'
        elif int(var_pos) <= 27759:
            region = 'ORF7A'
        elif int(var_pos) <= 27887:
            region = 'ORF7B'
        elif int(var_pos) <= 28259:
            region = 'ORF8'
        elif int(var_pos) <= 29533:
            region = 'N'
        elif int(var_pos) <= 29674:
            region = 'ORF10'
        else:
            region = '3UTR'
        if 'pos' in del_dict and int(line[1]) in del_dict['pos']:
            item_index = del_dict['pos'].index(int(line[1]))
            out_list.append([prefix,var_pos,region,int(var_depth)+int(del_dict['del_depth'][item_index]),A_depth,A_plus,A_minus,C_depth,C_plus,C_minus,G_depth,G_plus,G_minus,T_depth,T_plus,T_minus,del_dict['del_depth'][item_index],del_dict['del_depth_plus'][item_index],del_dict['del_depth_minus'][item_index],del_dict['deletion'][item_index]])
        elif len(line) > 10 and '+' in indel_seq:
            out_list.append([prefix,var_pos,region,int(var_depth)+int(indel_depth),A_depth,A_plus,A_minus,C_depth,C_plus,C_minus,G_depth,G_plus,G_minus,T_depth,T_plus,T_minus,indel_depth,indel_plus,indel_minus,indel_seq])
        else:
            out_list.append([prefix,var_pos,region,var_depth,A_depth,A_plus,A_minus,C_depth,C_plus,C_minus,G_depth,G_plus,G_minus,T_depth,T_plus,T_minus,0,0,0,'Na'])
    output_csv_writer.writerows(out_list)

with open(bam_rc_file+'.fmt.tsv','r') as bam_readcount_formated:
    df = pd.read_csv(bam_readcount_formated,sep='\t',header=0)
    conditions = [
        #A/T
        (df['A_DEPTH'] >= depth_value) & (df['T_DEPTH'] > depth_value) & (df['A_DEPTH']/df['DEPTH'] >= per_limit) & (df['T_DEPTH']/df['DEPTH'] >= per_limit) & (df['C_DEPTH']/df['DEPTH'] < per_limit) & (df['G_DEPTH']/df['DEPTH'] < per_limit) & ((df['A_PLUS']/df['A_MINUS'] > 0.05) & (df['INDEL_DEPTH']/df['DEPTH'] < per_limit) & (df['A_MINUS']/df['A_PLUS'] > 0.05)) & ((df['T_PLUS']/df['T_MINUS'] > 0.05) & (df['T_MINUS']/df['T_PLUS'] > 0.05)),
        #A/C
        (df['A_DEPTH'] >= depth_value) & (df['C_DEPTH'] > depth_value) & (df['A_DEPTH']/df['DEPTH'] >= per_limit) & (df['C_DEPTH']/df['DEPTH'] >= per_limit) & (df['T_DEPTH']/df['DEPTH'] < per_limit) & (df['G_DEPTH']/df['DEPTH'] < per_limit) & (df['INDEL_DEPTH']/df['DEPTH'] < per_limit) & ((df['A_PLUS']/df['A_MINUS'] > 0.05) & (df['A_MINUS']/df['A_PLUS'] > 0.05)) & ((df['C_PLUS']/df['C_MINUS'] > 0.05) & (df['C_MINUS']/df['C_PLUS'] > 0.05)),
        #A/G
        (df['A_DEPTH'] >= depth_value) & (df['G_DEPTH'] > depth_value) & (df['A_DEPTH']/df['DEPTH'] >= per_limit) & (df['G_DEPTH']/df['DEPTH'] >= per_limit) & (df['C_DEPTH']/df['DEPTH'] < per_limit) & (df['T_DEPTH']/df['DEPTH'] < per_limit) & (df['INDEL_DEPTH']/df['DEPTH'] < per_limit) & ((df['A_PLUS']/df['A_MINUS'] > 0.05) & (df['A_MINUS']/df['A_PLUS'] > 0.05)) & ((df['G_PLUS']/df['G_MINUS'] > 0.05) & (df['G_MINUS']/df['G_PLUS'] > 0.05)),
        #A/INDEL
        (df['A_DEPTH'] >= depth_value) & (df['INDEL_DEPTH'] > depth_value) & (df['A_DEPTH']/df['DEPTH'] >= per_limit) & (df['INDEL_DEPTH']/df['DEPTH'] >= per_limit) & ((df['A_PLUS']/df['A_MINUS'] > 0.05) & (df['A_MINUS']/df['A_PLUS'] > 0.05)) & ((df['INDEL_PLUS']/df['INDEL_MINUS'] > 0.05) & (df['INDEL_MINUS']/df['INDEL_PLUS'] > 0.05)),
        #T/C
        (df['T_DEPTH'] >= depth_value) & (df['C_DEPTH'] > depth_value) & (df['T_DEPTH']/df['DEPTH'] >= per_limit) & (df['C_DEPTH']/df['DEPTH'] >= per_limit) & (df['A_DEPTH']/df['DEPTH'] < per_limit) & (df['G_DEPTH']/df['DEPTH'] < per_limit) & (df['INDEL_DEPTH']/df['DEPTH'] < per_limit) & ((df['T_PLUS']/df['T_MINUS'] > 0.05) & (df['T_MINUS']/df['T_PLUS'] > 0.05)) & ((df['C_PLUS']/df['C_MINUS'] > 0.05) & (df['C_MINUS']/df['C_PLUS'] > 0.05)),
        #T/G
        (df['T_DEPTH'] >= depth_value) & (df['G_DEPTH'] > depth_value) & (df['T_DEPTH']/df['DEPTH'] >= per_limit) & (df['G_DEPTH']/df['DEPTH'] >= per_limit) & (df['A_DEPTH']/df['DEPTH'] < per_limit) & (df['C_DEPTH']/df['DEPTH'] < per_limit) & (df['INDEL_DEPTH']/df['DEPTH'] < per_limit) & ((df['T_PLUS']/df['T_MINUS'] > 0.05) & (df['T_MINUS']/df['T_PLUS'] > 0.05)) & ((df['G_PLUS']/df['G_MINUS'] > 0.05) & (df['G_MINUS']/df['G_PLUS'] > 0.05)),
        #T/INDEL
        (df['T_DEPTH'] >= depth_value) & (df['INDEL_DEPTH'] > depth_value) & (df['T_DEPTH']/df['DEPTH'] >= per_limit) & (df['INDEL_DEPTH']/df['DEPTH'] >= per_limit) & ((df['T_PLUS']/df['T_MINUS'] > 0.05) & (df['T_MINUS']/df['T_PLUS'] > 0.05)) & ((df['INDEL_PLUS']/df['INDEL_MINUS'] > 0.05) & (df['INDEL_MINUS']/df['INDEL_PLUS'] > 0.05)),
        #C/G
        (df['C_DEPTH'] >= depth_value) & (df['G_DEPTH'] > depth_value) & (df['C_DEPTH']/df['DEPTH'] >= per_limit) & (df['G_DEPTH']/df['DEPTH'] >= per_limit) & (df['A_DEPTH']/df['DEPTH'] < per_limit) & (df['T_DEPTH']/df['DEPTH'] < per_limit) & (df['INDEL_DEPTH']/df['DEPTH'] < per_limit) & ((df['C_PLUS']/df['C_MINUS'] > 0.05) & (df['C_MINUS']/df['C_PLUS'] > 0.05)) & ((df['G_PLUS']/df['G_MINUS'] > 0.05) & (df['G_MINUS']/df['G_PLUS'] > 0.05)),
        #C/INDEL
        (df['C_DEPTH'] >= depth_value) & (df['INDEL_DEPTH'] > depth_value) & (df['C_DEPTH']/df['DEPTH'] >= per_limit) & (df['INDEL_DEPTH']/df['DEPTH'] >= per_limit) & ((df['C_PLUS']/df['C_MINUS'] > 0.05) & (df['C_MINUS']/df['C_PLUS'] > 0.05)) & ((df['INDEL_PLUS']/df['INDEL_MINUS'] > 0.05) & (df['INDEL_MINUS']/df['INDEL_PLUS'] > 0.05)),
        #G/INDEL
        (df['G_DEPTH'] >= depth_value) & (df['INDEL_DEPTH'] > depth_value) & (df['G_DEPTH']/df['DEPTH'] >= per_limit) & (df['INDEL_DEPTH']/df['DEPTH'] >= per_limit) & ((df['G_PLUS']/df['G_MINUS'] > 0.05) & (df['G_MINUS']/df['G_PLUS'] > 0.05)) & ((df['INDEL_PLUS']/df['INDEL_MINUS'] > 0.05) & (df['INDEL_MINUS']/df['INDEL_PLUS'] > 0.05)),
        #A/T/C
        (df['A_DEPTH'] >= depth_value) & (df['T_DEPTH'] > depth_value) & (df['C_DEPTH'] > depth_value) & (df['A_DEPTH']/df['DEPTH'] >= per_limit) & (df['T_DEPTH']/df['DEPTH'] >= per_limit) & (df['C_DEPTH']/df['DEPTH'] >= per_limit) & (df['G_DEPTH']/df['DEPTH'] < per_limit) & ((df['A_PLUS']/df['A_MINUS'] > 0.05) & (df['A_MINUS']/df['A_PLUS'] > 0.05)) & ((df['T_PLUS']/df['T_MINUS'] > 0.05) & (df['T_MINUS']/df['T_PLUS'] > 0.05)) & ((df['C_PLUS']/df['C_MINUS'] > 0.05) & (df['C_MINUS']/df['C_PLUS'] > 0.05)),
        #A/T/G
        (df['A_DEPTH'] >= depth_value) & (df['T_DEPTH'] > depth_value) & (df['G_DEPTH'] > depth_value) & (df['A_DEPTH']/df['DEPTH'] >= per_limit) & (df['T_DEPTH']/df['DEPTH'] >= per_limit) & (df['G_DEPTH']/df['DEPTH'] >= per_limit) & (df['C_DEPTH']/df['DEPTH'] < per_limit) & ((df['A_PLUS']/df['A_MINUS'] > 0.05) & (df['A_MINUS']/df['A_PLUS'] > 0.05)) & ((df['T_PLUS']/df['T_MINUS'] > 0.05) & (df['T_MINUS']/df['T_PLUS'] > 0.05)) & ((df['G_PLUS']/df['G_MINUS'] > 0.05) & (df['G_MINUS']/df['G_PLUS'] > 0.05)),
        #A/G/C
        (df['A_DEPTH'] >= depth_value) & (df['C_DEPTH'] > depth_value) & (df['G_DEPTH'] > depth_value) & (df['A_DEPTH']/df['DEPTH'] >= per_limit) & (df['C_DEPTH']/df['DEPTH'] >= per_limit) & (df['G_DEPTH']/df['DEPTH'] >= per_limit) & (df['T_DEPTH']/df['DEPTH'] < per_limit) & ((df['A_PLUS']/df['A_MINUS'] > 0.05) & (df['A_MINUS']/df['A_PLUS'] > 0.05)) & ((df['G_PLUS']/df['G_MINUS'] > 0.05) & (df['G_MINUS']/df
        ['G_PLUS'] > 0.05)) & ((df['C_PLUS']/df['C_MINUS'] > 0.05) & (df['C_MINUS']/df['C_PLUS'] > 0.05)),
        #T/G/C
        (df['T_DEPTH'] >= depth_value) & (df['G_DEPTH'] > depth_value) & (df['C_DEPTH'] > depth_value) & (df['T_DEPTH']/df['DEPTH'] >= per_limit) & (df['G_DEPTH']/df['DEPTH'] >= per_limit) & (df['C_DEPTH']/df['DEPTH'] >= per_limit) & (df['A_DEPTH']/df['DEPTH'] < per_limit) & ((df['T_PLUS']/df['T_MINUS'] > 0.05) & (df['T_MINUS']/df['T_PLUS'] > 0.05)) & ((df['G_PLUS']/df['G_MINUS'] > 0.05) & (df['G_MINUS']/df['G_PLUS'] > 0.05)) & ((df['C_PLUS']/df['C_MINUS'] > 0.05) & (df['C_MINUS']/df['C_PLUS'] > 0.05)),
        #A/T/G/C        
        (df['A_DEPTH'] >= depth_value) & (df['T_DEPTH'] > depth_value) & (df['C_DEPTH'] > depth_value) & (df['G_DEPTH'] > depth_value) & (df['A_DEPTH']/df['DEPTH'] >= per_limit) & (df['T_DEPTH']/df['DEPTH'] >= per_limit) & (df['C_DEPTH']/df['DEPTH'] >= per_limit) & (df['G_DEPTH']/df['DEPTH'] >= per_limit) & ((df['A_PLUS']/df['A_MINUS'] > 0.05) & (df['A_MINUS']/df['A_PLUS'] > 0.05)) & ((df['T_PLUS']/df['T_MINUS'] > 0.05) & (df['T_MINUS']/df['T_PLUS'] > 0.05)) & ((df['C_PLUS']/df['C_MINUS'] > 0.05) & (df['C_MINUS']/df['C_PLUS'] > 0.05)) & ((df['G_PLUS']/df['G_MINUS'] > 0.05) & (df['G_MINUS']/df['G_PLUS'] > 0.05))]
    minors = ['A/T','A/C','A/G','A/INDEL','T/C','T/G','T/INDEL','C/G','C/INDEL','G/INDEL','A/T/C','A/T/G','A/G/C','T/G/C','A/T/C/G']
    df['PUTATIVE_MINOR'] = np.select(conditions,minors,default='False')
    df = df.loc[df['PUTATIVE_MINOR'] != "False"]
    df.to_csv(bam_rc_file+'.fmt.minors.tsv',sep='\t',index=False)
os.remove(bam_rc_file+'.fmt.tsv')