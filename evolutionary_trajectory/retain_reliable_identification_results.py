#!/usr/bin/env python
# -*- coding: utf-8 -*-
#writed by Lhx on 4/15/2022
import argparse
parser = argparse.ArgumentParser(description="Retain reliable identification results!")
parser.add_argument("-CLOUD",required=True,help="CLOUD final results")
parser.add_argument("-CDROM",required=True,help="CDROM results")
parser.add_argument("-binary",required=True,help="binary results")
parser.add_argument("-o",required=True,help="output prefix")
args = parser.parse_args()

outfile = open(args.o+'.classifications','w')
d={}
print('Parent','Child','Ancestor','Classification',sep='\t',file=outfile)
with open(args.CLOUD) as f:
    for lines in f.readlines()[1:]:
        line = lines.strip().split()
        d[line[0]+'\t'+line[1]+'\t'+line[2]] = [line[3]]
with open(args.CDROM) as f:
    for lines in f.readlines()[1:]:
        line = lines.strip().split()
        d[line[0]+'\t'+line[1]+'\t'+line[2]].append(line[6])
with open(args.binary) as f:
    for lines in f.readlines()[1:]:
        line = lines.strip().split()
        d[line[0]+'\t'+line[1]+'\t'+line[2]].append(line[6])
l=[]
for k,v in d.items():
    if v[0] == v[1] and v[0] == v[2]:
        print(k,v[0],sep='\t',file=outfile)
        l.append(v[0])
    elif v[0] == v[1] and v[0] != v[2]:
        print(k,v[0],sep='\t',file=outfile)
        l.append(v[0])
    elif v[0] == v[2] and v[0] != v[1]:
        print(k,v[0],sep='\t',file=outfile)
        l.append(v[0])
    elif v[1] == v[2] and v[0] != v[1]:
        print(k,v[1],sep='\t',file=outfile)
        l.append(v[1])
    elif v[0] != v[1] and v[0] != v[2] and v[1] != v[2]:
        print(k,v[0],sep='\t',file=outfile)
        l.append(v[0])
outfile.close()
count_file = open(args.o+'.counts','w')
print('Conservation', l.count('Conservation'), sep='\t',file=count_file)
print('Neofunctionalization(Parent)', l.count('Neofunctionalization(Parent)'), sep='\t', file=count_file)
print('Neofunctionalization(Child)', l.count('Neofunctionalization(Child)'), sep='\t', file=count_file)
print('Subfunctionalization', l.count('Subfunctionalization'), sep='\t', file=count_file)
print('Specialization', l.count('Specialization'), sep='\t',file=count_file)
count_file.close()
