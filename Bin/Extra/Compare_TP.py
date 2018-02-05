##################
##################  MERGE POST-PROCCESING
##################
#
#
# The SNVs for which its VSR is classified as a TP region for all three individuals are identified.
#
# COMMAND LINE.
# python Compare_TP.py -f HG1_TP.tab -m HG2_TP.tab -c HG3_TP.tab -o Novo_TP.tab
#
# PARAMETERS
#	f 		File with the TP de novo SNVs for all chromosomes of the father
#	m 		File with the TP de novo SNVs for all chromosomes of the mother 
#	c		File with the TP de novo SNVs for all chromosomes of the child
#	out 		Output file
#
# OUTPUT
# Same columns as Novo.out
################

import sys, os, getopt
import argparse


def obtain_data(file, array, data_dict):
    with open(file) as handle:
        lines = handle.readlines()

    for line in lines:
        line.rstrip()
        info = line.split("\t")
        chr = info[0]
        pos = info[4]
        data = chr + "-" + pos
        array.append(data)
        data_dict[data] = line
    return array


parser = argparse.ArgumentParser()
parser.add_argument('-f', '--fatherTP ', help="File with the whole-genome TP de novo SNVs for the father", dest='father', required=True, metavar="father_TP")
parser.add_argument("-m", "--motherTP", help="File with the whole-genome TP de novo SNVs for the mother  ", dest='mother',required=True, metavar="mother_TP")
parser.add_argument("-c", "--childTP", help="File with the whole-genome TP de novo SNVs for the child", dest='child',required=True, metavar="child_TP")
parser.add_argument("-o", "--output", help="Output file. TP in all three individuals",dest='out',required=True, metavar="all_TP")
parser.add_argument("-p", "--outputALL", help="More information about TP in all three individuals", dest='out_all',required=True,metavar="all_TP_detail")

args = parser.parse_args()

father_file = args.father
mother_file = args.mother
child_file = args.child
out_file = args.out
out_all_file = args.out_all

        
out = open(out_file, 'w')
out_all= open(out_all_file, 'w')

data_dict = {}

father_tp = []
mother_tp = []
child_tp = []
##### READ FATHER TP FILE   
father_tp = obtain_data(father_file, father_tp, data_dict)

##### READ MOTHER TP FILE   
mother_tp = obtain_data(mother_file, mother_tp, data_dict)

##### READ CHILD TP FILE   
child_tp = obtain_data(child_file, child_tp, data_dict)


for tp in child_tp:
    if tp in father_tp and tp in mother_tp:
        out.write(tp + '\n')
        out_all.write(data_dict[tp])

