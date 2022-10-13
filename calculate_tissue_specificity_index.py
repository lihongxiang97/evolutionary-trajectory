#! /usr/bin/env python
#Writed by LiHongxiang on 3/24/2022
#根据不同组织的表达量，计算基因组织特异性系数
import argparse
import math
parser = argparse.ArgumentParser(description="Calculate the tissue specificity index!")
parser.add_argument("-T",required=True,help="Gene exp file:format:gene_id\ttissue1_exp\ttissue2_exp...\nPbr00001\t3\t5...")
parser.add_argument("-N",required=True,help="Number of tissues")
parser.add_argument("-o",required=True,help="outfile,format:gene_id\tindex_value)
args = parser.parse_args()
exp_file = args.T
n = args.N
outfile = args.o
#读文件
def read_f(file):
    f = open(file).readlines()
    return f
exp = read_f(exp_file)

#计算基因组织特异性系数
d={}
for lines in exp:
    line = lines.strip().split()
    max = 0
    for i in range(1,int(n)+1):
        if float(line[i]) >= max:
            max = float(line[i])
    index = 0
    for i in range(1,int(n)+1):
        index += (1-math.log(float(line[i])+1)/math.log(max+1))/(int(n)-1)
    d[line[0]] = index
out = open(outfile,'w')
for key,value in d.items():
    print(key,value,sep='\t',file=out)




