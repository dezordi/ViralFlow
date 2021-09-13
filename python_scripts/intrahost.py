#!/usr/bin/python3
# -*- coding: utf-8 -*-

__author__ = "Filipe Z. Dezordi"
__license__ = "GPL"
__version__ = "0.0.1"
__maintainer__ = "Filipe Z. Dezordi"
__email__ = "zimmer.filipe@gmail.com"
__date__ = "2021/08/06"
__username__ = "dezordi"

import argparse, csv, re, os
import pandas as pd
import numpy as np
from Bio import AlignIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO

parser = argparse.ArgumentParser(description = 'This scripts returns putative minor variants, based on nucleotide diversity of reads mapped against a reference genome')
parser.add_argument("-in","--bam_readcount",help="bam-readcount output file", required=True)
parser.add_argument("-lm","--per_limit",help="Percentage limit of Depth Base / Total Base to consider as putative minor variant, default = 0.05, which represents 5 percent", type=float,default= 0.05)
parser.add_argument("-dp","--depth_number",help="Minimum depth value to consider as minor variant, default = 100", type=int, default= 100)
parser.add_argument("-al","--alignment_file",help="Alignment file against wuhan reference genome, with keeplength parameter", required=True)
args = parser.parse_args()
#store arguments
bam_rc_file = args.bam_readcount
per_limit = args.per_limit
depth_value = args.depth_number
alignment_file = args.alignment_file
#create prefix
prefix = bam_rc_file
prefix = re.sub(".*\/","",prefix)
prefix = re.sub("\..*","",prefix)
prefix = re.sub(r'__','/',prefix)
prefix = re.sub(r'--','|',prefix)
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
        A_depth = int(re.sub(r".*:","",A_depth))
        A_plus = ":".join(A_depth_line.split(":", 6)[:6])
        A_plus = int(re.sub(r".*:","",A_plus))
        A_minus = ":".join(A_depth_line.split(":", 7)[:7])
        A_minus = int(re.sub(r".*:","",A_minus))
        C_depth_line = line[6].rstrip('\n')
        C_depth = ":".join(C_depth_line.split(":", 2)[:2])
        C_depth = int(re.sub(r".*:","",C_depth))
        C_plus = ":".join(C_depth_line.split(":", 6)[:6])
        C_plus = int(re.sub(r".*:","",C_plus))
        C_minus = ":".join(C_depth_line.split(":", 7)[:7])
        C_minus = int(re.sub(r".*:","",C_minus))
        G_depth_line = line[7].rstrip('\n')
        G_depth = ":".join(G_depth_line.split(":", 2)[:2])
        G_depth = int(re.sub(r".*:","",G_depth))
        G_plus = ":".join(G_depth_line.split(":", 6)[:6])
        G_plus = int(re.sub(r".*:","",G_plus))
        G_minus = ":".join(G_depth_line.split(":", 7)[:7])
        G_minus = int(re.sub(r".*:","",G_minus))
        T_depth_line = line[8].rstrip('\n')
        T_depth = ":".join(T_depth_line.split(":", 2)[:2])
        T_depth = int(re.sub(r".*:","",T_depth))
        T_plus = ":".join(T_depth_line.split(":", 6)[:6])
        T_plus = int(re.sub(r".*:","",T_plus))
        T_minus = ":".join(T_depth_line.split(":", 7)[:7])
        T_minus = int(re.sub(r".*:","",T_minus))
        major_depth_column = 0
        depth_dict = {A_depth:"A_depth",C_depth:"C_depth",G_depth:"G_depth",T_depth:"T_Depth"}
        major_depth = max(depth_dict)
        major_depth_plus = max([A_plus,C_plus,G_plus,T_plus])
        major_depth_minus = max([A_minus,C_minus,G_minus,T_minus])     
        if len(line) > 11: #bamreadcount output can have 11 or more columns, depending of number of putative indels
            for i in range(10,len(line)):
                indel_depth_line = line[i].rstrip('\n')
                indel_seq = re.sub(r":.*","",indel_depth_line)
                indel_depth = ":".join(indel_depth_line.split(":", 2)[:2])
                indel_depth = int(re.sub(r".*:","",indel_depth))
                if '+' in indel_seq:
                    major_depth = int(major_depth) - int(indel_depth)
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
            if '+' in indel_seq:
                if depth_dict.get(max(depth_dict)) == 'A_depth':
                    A_depth = int(A_depth) - int(indel_depth)
                    A_plus = int(A_plus) - int(indel_plus)
                    A_minus = int(A_minus) - int(indel_minus)
                elif depth_dict.get(max(depth_dict)) == 'C_depth':
                    C_depth = int(C_depth) - int(indel_depth)
                    C_plus = int(C_plus) - int(indel_plus)
                    C_minus = int(C_minus) - int(indel_minus)
                elif depth_dict.get(max(depth_dict)) == 'G_depth':
                    G_depth = int(G_depth) - int(indel_depth)
                    G_plus = int(G_plus) - int(indel_plus)
                    G_minus = int(G_minus) - int(indel_minus)
                elif depth_dict.get(max(depth_dict)) == 'T_depth':
                    T_depth = int(T_depth) - int(indel_depth)
                    T_plus = int(T_plus) - int(indel_plus)
                    T_minus = int(T_minus) - int(indel_minus)
            elif '-' in indel_seq:
                if int(indel_depth) < depth_value:
                    pass
                else:
                    del_length = len(indel_seq)-1
                    for i in range(0,del_length):
                        if i == 0:
                            var_depth = int(var_depth) - int(indel_depth)
                            del_dict['pos'].append(int(line[1])+i)
                            del_dict['deletion'].append(indel_seq)
                            del_dict['del_depth'].append(indel_depth)
                            del_dict['del_depth_plus'].append(indel_plus)
                            del_dict['del_depth_minus'].append(indel_minus)
                        else:
                            del_dict['pos'].append(int(line[1])+i)
                            del_dict['deletion'].append(indel_seq)
                            del_dict['del_depth'].append(indel_depth)
                            del_dict['del_depth_plus'].append(indel_plus)
                            del_dict['del_depth_minus'].append(indel_minus)
        elif len(line) == 11:
            indel_depth_line = line[10].rstrip('\n')
            indel_seq = re.sub(r":.*","",indel_depth_line)
            indel_depth = ":".join(indel_depth_line.split(":", 2)[:2])
            indel_depth = int(re.sub(r".*:","",indel_depth))
            indel_plus = ":".join(indel_depth_line.split(":", 6)[:6])
            indel_plus = int(re.sub(r".*:","",indel_plus))
            indel_minus = ":".join(indel_depth_line.split(":", 7)[:7])
            indel_minus = int(re.sub(r".*:","",indel_minus))
            if '+' in indel_seq:
                major_depth = int(major_depth) - int(indel_depth)
                if depth_dict.get(max(depth_dict)) == 'A_depth':
                    A_depth = int(A_depth) - int(indel_depth)
                    A_plus = int(A_plus) - int(indel_plus)
                    A_minus = int(A_minus) - int(indel_minus)
                elif depth_dict.get(max(depth_dict)) == 'C_depth':
                    C_depth = int(C_depth) - int(indel_depth)
                    C_plus = int(C_plus) - int(indel_plus)
                    C_minus = int(C_minus) - int(indel_minus)
                elif depth_dict.get(max(depth_dict)) == 'G_depth':
                    G_depth = int(G_depth) - int(indel_depth)
                    G_plus = int(G_plus) - int(indel_plus)
                    G_minus = int(G_minus) - int(indel_minus)
                elif depth_dict.get(max(depth_dict)) == 'T_depth':
                    T_depth = int(T_depth) - int(indel_depth)
                    T_plus = int(T_plus) - int(indel_plus)
                    T_minus = int(T_minus) - int(indel_minus)
            if '-' in indel_seq:
                if int(indel_depth) < depth_value:
                    pass
                else:
                    del_length = len(indel_seq)-1
                    for i in range(0,del_length):
                        if i == 0:
                            var_depth = int(var_depth) - int(indel_depth)
                            del_dict['pos'].append(int(line[1])+i)
                            del_dict['deletion'].append(indel_seq)
                            del_dict['del_depth'].append(indel_depth)
                            del_dict['del_depth_plus'].append(indel_plus)
                            del_dict['del_depth_minus'].append(indel_minus)
                        else:
                            del_dict['pos'].append(int(line[1])+i)
                            del_dict['deletion'].append(indel_seq)
                            del_dict['del_depth'].append(indel_depth)
                            del_dict['del_depth_plus'].append(indel_plus)
                            del_dict['del_depth_minus'].append(indel_minus)
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
            out_list.append([prefix,var_pos,region,var_depth,A_depth,A_plus,A_minus,C_depth,C_plus,C_minus,G_depth,G_plus,G_minus,T_depth,T_plus,T_minus,indel_depth,indel_plus,indel_minus,indel_seq])
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

#Formating intrahost variants

output_tsv = open(bam_rc_file+'.intrahost.tsv','w')
output_tsv_writer = csv.writer(output_tsv, delimiter='\t')
output_tsv_writer.writerow(["GENOME","POS","REGION","DEPTH","A_DEPTH","A_PLUS","A_MINUS","C_DEPTH","C_PLUS","C_MINUS","G_DEPTH","G_PLUS","G_MINUS","T_DEPTH","T_PLUS","T_MINUS","INDEL_DEPTH","INDEL_PLUS","INDEL_MINUS","INDEL","PUTATIVE_MINOR","MAJOR","MINOR","MAJOR_DEPTH","MINOR_DEPTH","FREQUENCY"])
out_list_lists = []
output_short_tsv = open(bam_rc_file+'.intrahost.short.tsv','w')
output_short_writer = csv.writer(output_short_tsv, delimiter='\t')
output_short_writer.writerow(["GENOME","POS","REGION","DEPTH","MAJOR","MINOR","MAJOR_DEPTH","MINOR_DEPTH","FREQUENCY"])
out_short_list_lists = []
with open(bam_rc_file+'.fmt.minors.tsv','r') as minor_tsv:
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
        if 'INDEL' not in line[20]:
            line[19] = 'Na'
        if 'GENOME' not in line[0]:
            try: #For intrahost loci with 3 or more alleles, don't check frequency
                frequency = "{:.2f}".format(int(minor_depth)/int(line[3].rstrip('\n')))
            except:
                frequency = 'CURATION'
            out_list_lists.append([line[0],line[1],line[2],line[3],line[4],line[5],line[6],line[7],line[8],line[9],line[10],line[11],line[12],line[13],line[14],line[15],line[16],line[17],line[18],line[19],line[20],major,minor,major_depth,minor_depth,frequency])
            out_short_list_lists.append([line[0],line[1],line[2],line[3],major,minor,major_depth,minor_depth,frequency])
output_tsv_writer.writerows(out_list_lists)
output_short_writer.writerows(out_short_list_lists)
output_tsv.close()
#generate minor consensus
output_handle_min = open(alignment_file +'.minor.fa', "w")
alignment = AlignIO.read(open(alignment_file),'fasta')

for record in alignment:
    new_seq_lst = list(record.seq)
    new_seq = ''
    record_lst = []
    with open(bam_rc_file+'.intrahost.tsv','r') as minor_variants_file:
        minor_reader = csv.reader(minor_variants_file,delimiter='\t')
        for line in minor_reader:
            line[0] = re.sub(r'__','/',line[0])
            line[0] = re.sub(r'--','|',line[0])
            if record.id == line[0]:
                pos = int(line[1].rstrip('\n'))
                minor = line[22].rstrip('\n')
                new_seq_lst[pos-1] = minor
    for letter in new_seq_lst:
        new_seq += letter
    if record.id != 'NC_045512.2':
        record = SeqRecord(Seq(new_seq), id=record.id+'_minor', description = '')
        record_lst.append(record)
        SeqIO.write(record_lst,output_handle_min,'fasta')

os.remove(bam_rc_file+'.fmt.minors.tsv')
os.remove(bam_rc_file+'.fmt.tsv')