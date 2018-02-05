#!/bin/python


################
################ OBTAIN SIGNATURE CS'S SEQUENCE
################  
#
# This script will obtain the sequences for the PreCS and PostCS for every VSR.
# All the VSR file must be located on the same directory and they must have the same extension
# The header of the multi-FASTA file that is generated will change for the CHILD and for the PARENTS in the family-based framework
# 
#             COMMAND LINE.
#     python Cut_SignatureCSs.py --VSR=SR/ --REFDIR=REFDIR/ --sufixREF=sufixREF --sufixVSR=sufixVSR --prefixVSR=prefixVSR --kSIZE=kmer_size --output=SignCS.fa
#     python Cut_SignatureCSs.py -v SR/ -r REFDIR/ -x sufixREF -y sufixVSR -z prefixVSR -k kmer_size -o SignCS.fa
#
# PARAMETERS.
#	Script 		Cut_SignatureCSs.py
#	VSR 		VSR directory path
#	REFDIR 		RG directory path
#	sufixREF 	Sufix of the fasta reference files
#	prefixVSR 	Prefix of VSR files (if the VSR files have no prefix, set to NA)
#	sufixVSR 	Sufix of VSR files
#	kmer_size	CS size (SET to 30)
#	SignCS 		Output file
#
# The name of the files must follow the next rules. Example:
# RG chr1 file name: chr1.{sufixREF}
# VSR chr1 file name: {prefixSR}chr1{sufixSR}
#
# OUTPUT.
# A multi-fasta file with the sequences of all Signature CSs
#
################
################ 
################

import sys 
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-v', '--vsr ', help="VSR directory path", dest='vsr', required=True, metavar="SR/")
parser.add_argument("-r", "--REFDIR", help="RG directory path", dest='REFDIR',required=True, metavar="REFDIR/")
parser.add_argument("-x", "--sufixREF", help="Sufix of the fasta reference files",dest='sufixREF',required=True, metavar="sufixREF")
parser.add_argument("-y", "--sufixVSR", help="Sufix of VSR files", dest='sufixVSR',required=True,metavar="sufixVSR")
parser.add_argument("-z", "--prefixVSR", help="Prefix of VSR files (if the VSR files have no prefix, set to NA)", dest='prefixVSR',required=True,metavar="prefixVSR")
parser.add_argument("-k", "--kSIZE", help="CS size (SET to 30)", dest='kSIZE',required=True,metavar="kmer_size")
parser.add_argument("-o", "--output", help="Output file", dest='output',required=True,metavar="SignCS")

args = parser.parse_args()

dir_var = args.vsr
dir_fna = args.REFDIR
sufix_fna = args.sufixREF
sufix_var = args.sufixVSR
prefix_var = args.prefixVSR
k = args.kSIZE
out = args.output


#(script, dir_var, dir_fna, sufix_fna, prefix_var, sufix_var, k, out)= sys.argv


#### THIS WILL STORE ALL VSR FILES
filenames = os.listdir(dir_var)
var_filenames = [s for s in filenames if s.endswith(sufix_var)]

print (var_filenames)

final = []

#### PROCESS EVERY VSR FILES
for var in var_filenames:
    chr = var.replace(prefix_var, "")
    chr = chr.replace(sufix_var, "")
    fna = chr + "." + sufix_fna	
    fna = dir_fna + fna	
    var = dir_var + var
	
	#### WILL READ FASTA FILE
    handle = open(fna, "rU")
    fasta = SeqIO.read(handle, "fasta")
    handle.close()

    #### WILL LOAD VSR FILE INTO MEMORY
    var_file = open(var, 'r').read()
    var_lines = var_file.split('\n')
    if (var_lines[-1] == ""):
        var_lines = var_lines[:-1]


    #### PROCESSING EVERY VSR
    for line in var_lines:
        info = line.split('\t')
        if info[0] == '':
            continue
                    
        #### LOAD PRE-CS AND POST-CS START POSITIONS
        cs_prev = int(info[0])
        cs_after = int(info[8])

        #### CALCULATE PRE-CS AND POST-CS END POSITIONS
        st_prev = cs_prev - 1
        end_prev = st_prev + int(k)
        st_after = cs_after - 1
        end_after = st_after + int(k)
	

	#### CUT PRE-CS AND POST-CS SEQUENCE FROM REFERENCE 
	 #### APPEND SEQUENCES TO MILTI-FASTA OBJECT
        sub = fasta.seq[st_prev:end_prev]
        record = SeqRecord(sub, '%s:%i:%i:%i' % (chr, cs_prev, cs_prev, cs_after), '', '')
        final.append(record)
        sub = fasta.seq[st_after:end_after]
        record = SeqRecord(sub, '%s:%i:%i:%i' % (chr, cs_after, cs_prev, cs_after), '', '')
        final.append(record)

#print final

output_handle = open(out, "w")
SeqIO.write(final, output_handle, "fasta")
output_handle.close()

