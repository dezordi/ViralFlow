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
from Bio import AlignIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO

parser = argparse.ArgumentParser(description = 'This script replace major to minor variants')
parser.add_argument("-in","--input_file",help="Alignment file against wuhan reference genome", required=True)
parser.add_argument("-mv","--minor_variants",help="TSV minor formated file", required=True)
args = parser.parse_args()
input_file = args.input_file
minor_variants = args.minor_variants
output_handle_min = open(input_file+'.minor.fa', "w")
alignment = AlignIO.read(open(input_file),'fasta')

for record in alignment:
    new_seq_lst = list(record.seq)
    new_seq = ''
    record_lst = []
    with open(minor_variants,'r') as minor_variants_file:
        minor_reader = csv.reader(minor_variants_file,delimiter='\t')
        for line in minor_reader:
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