#!/bin/python

##################
##################
## OBTAIN SIGNATURE CS'S SEQUENCE (PARENTS)
##################
##################
#
#		 python SignatureCS_Parent.py -c SignatureCSs_child.fa -o SignatureCSs_parent.fa
# 	python SignatureCS_Parent.py --SignatureCSs SignatureCSs_child.fa --OUTPUT SignatureCSs_parent.fa
#
# PARAMETERS.
#	Script 				SignatureCS_Parent.py
#	SignatureCSs, c		Multi-fasta file with all Signature CSs for the child
#	OUTPUT, o	 		Output file (multi-fasta file), the same SginatureCS with IDs for the parents
#
#
# OUTPUT.
# A multi-fasta file with the sequences of all retrieved CSs sequence.
# The identifier if every fasta sequence contains:
#   -CS ID for the sequence retrieved
#   -PrevCS for the child VSR
#   -Post CS for the child VSR
# The last two positions will ve repeated twice
#
#
# NOTES.
# This step requires only one process
##################
##################


import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument('-c', '--SignatureCSs ', help="Multi-fasta file with all Signature CSs for the child", dest='CSs', required=True, metavar="SignatureCSs_child.fa")
parser.add_argument("-o", "--OUTPUT", help="Output file (multi-fasta file), the same SginatureCS with IDs for the parents", dest='output',required=True, metavar="SignatureCSs_parent.fa/")

args = parser.parse_args()

cs_child = args.CSs
out = args.output



#import sys 

#(script, cs_child, out)= sys.argv

out_file = open(out, 'w')
cs_file = open(cs_child)


for line in cs_file:
    line = line.rstrip('\n')
    if '>' in line:
        info = line.split(':')
        out_file.write(line + ":" +  info[2] + ":" + info[3] + "\n")
    else:
        out_file.write(line + "\n")

