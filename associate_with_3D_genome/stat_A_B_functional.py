#! /usr/bin/env python
# -*- coding: utf-8 -*-
#Writed by LiHongxiang on 2/28/2022
import argparse
parser = argparse.ArgumentParser(description="Statistics A/B compartment and their transformations!")
parser.add_argument("-C",required=True,help="functionalization_compartment")
parser.add_argument("-o",required=True,help="outfile name")
args = parser.parse_args()
out = args.o
c = args.C
outfile = open(out,'w')
compartment = open(c).readlines()

def stat(file,type):
    A_A = 0
    B_B = 0
    A_B = 0
    B_A = 0
    l=[]
    for lines in file:
        line = lines.strip().split()
        if line[4] == type:
            if line[1] == line[3] == 'A':
                A_A += 1
            if line[1] == line[3] == 'B':
                B_B += 1
            if line[1] == 'A' and line[3] == 'B':
                A_B += 1
            if line[1] == 'B' and line[3] == 'A':
                B_A += 1
    l.append(['A-A', str(A_A), type])
    l.append(['B-B', str(B_B), type])
    l.append(['A-B', str(A_B), type])
    l.append(['B-A', str(B_A), type])
    return l
con = stat(compartment,'Conservation')
spe = stat(compartment,'Specialization')
neo_p = stat(compartment,'Neofunctionalization(Parent)')
neo_c = stat(compartment,'Neofunctionalization(Child)')
sub = stat(compartment,'Subfunctionalization')

l = con + spe + neo_p + neo_c + sub
for i in l:
    print('\t'.join(i),file=outfile)



