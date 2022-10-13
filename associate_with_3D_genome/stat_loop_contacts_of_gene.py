#! /usr/bin/env python
# -*- coding: utf-8 -*-
#Writed by LiHongxiang on 3/19/2022
import argparse
parser = argparse.ArgumentParser(description="Statistics of GGL and DGL contacts number associate gene!")
parser.add_argument("-L",required=True,help="gene_gene loop or dOCR-gene loop file,out by annotation_loop.py")
parser.add_argument("-D",required=True,help="gene file,format:gene_id")
parser.add_argument("-G",help="gff file out by prepare.pl")
parser.add_argument("-o",required=True,help="outfile name")
parser.add_argument("-r",required=True,help="resolution")
args = parser.parse_args()
r = int(args.r)

out = open(args.o,'w')
loop = open(args.L).readlines()
dup = open(args.D).readlines()
contact={}
location={}
gene_list = []
for lines in dup:
    line = lines.strip().split()
    gene_list.append(line[0])
gene_list_uniq = list(set(gene_list))
for i in gene_list_uniq:
    for rows in loop:
        if i in rows:
            row = rows.strip().split()
            #制作基因名对应contacts数的字典
            if i in contact:
                contact[i] += int(row[4])
            else:
                contact[i] = int(row[4])
            #制作基因名对应位置的字典
            if i in location:
                location[i].append((row[0],row[1]))
                location[i].append((row[2],row[3]))
            else:
                location[i] = [(row[0],row[1]),(row[2],row[3])]
    try:
        location[i] = list(set(location[i]))
    except KeyError:
        pass
#基因位置的字典
d_gff={}
gff = open(args.G).readlines()
for lines in gff:
    line = lines.strip().split()
    if '-' in line[0]:
        chr = line[0].split('-')[1]
    else:
        chr = line[0]
    d_gff[line[1]] = [chr,line[2],line[3]]

#计算基因和loop重叠的长度
gene_length = {}
loop_length = {}
for id,sets in location.items():
    for l in sets:
        if l[0] == d_gff[id][0]:
            if id not in gene_length:
                if int(d_gff[id][1]) <= int(l[1]) - r/2 <= int(l[1]) + r/2 <= int(d_gff[id][2]):
                    gene_length[id] = r/1000
                    loop_length[id] = r/1000
                if int(l[1]) - r/2 <= int(d_gff[id][1]) <= int(d_gff[id][2]) <= int(l[1]) + r/2:
                    gene_length[id] = (int(d_gff[id][2]) - int(d_gff[id][1])) / 1000
                    loop_length[id] = r/1000
                if int(d_gff[id][1]) <= int(l[1]) - r/2 <= int(d_gff[id][2]) <= int(l[1]) + r/2:
                    gene_length[id] = (int(d_gff[id][2]) - int(l[1]) + r/2) / 1000
                    loop_length[id] = r/1000
                if int(l[1]) - r/2 <= int(d_gff[id][1]) <= int(l[1]) + r/2 <= int(d_gff[id][2]):
                    gene_length[id] = (int(l[1]) + r/2 - int(d_gff[id][1])) / 1000
                    loop_length[id] = r/1000
            else:
                if int(d_gff[id][1]) <= int(l[1]) - r/2 <= int(l[1]) + r/2 <= int(d_gff[id][2]):
                    gene_length[id] += r/1000
                    loop_length[id] += r/1000
                if int(l[1]) - r/2 <= int(d_gff[id][1]) <= int(d_gff[id][2]) <= int(l[1]) + r/2:
                    gene_length[id] += (int(d_gff[id][2]) - int(d_gff[id][1])) / 1000
                    loop_length[id] += r/1000
                if int(d_gff[id][1]) <= int(l[1]) - r/2 <= int(d_gff[id][2]) <= int(l[1]) + r/2:
                    gene_length[id] += (int(d_gff[id][2]) - int(l[1]) + r/2) / 1000
                    loop_length[id] += r/1000
                if int(l[1]) - r/2 <= int(d_gff[id][1]) <= int(l[1]) + r/2 <= int(d_gff[id][2]):
                    gene_length[id] += (int(l[1]) + r/2 - int(d_gff[id][1])) / 1000
                    loop_length[id] += r/1000

#基因占loop的比例乘loop上的contacts数，得到基因所占的contacts数
for key,value in contact.items():
    number = int(value)*(float(gene_length[key]) / int(loop_length[key]))
    contact[key] = str(number)

for i in gene_list:
    if i in contact:
        print(i,contact[i],sep='\t',file=out)
    else:
        print(i, '0', sep='\t', file=out)
