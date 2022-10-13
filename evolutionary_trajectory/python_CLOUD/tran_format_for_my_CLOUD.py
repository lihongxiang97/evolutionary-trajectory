#! /usr/bin/env python
# -*- coding: utf-8 -*-
#Writed by LiHongxiang on 4/8/2022
import argparse
parser = argparse.ArgumentParser(description="tran format for CLOUD!")
parser.add_argument("-T",required=True,help="divergencetime file")
parser.add_argument("-A",required=True,help="ancestral species TPM")
parser.add_argument("-B",required=True,help="target species TPM")
parser.add_argument("-o",required=True,help="outfile name")
parser.add_argument("-Singledata", type= bool, help="If singledata,set 1")
parser.add_argument("-Duplicatedata", type=bool, help="If duplicatedata,set 1")
args = parser.parse_args()
out = args.o
outfile = open(out,'w')
single = args.Singledata
duplicate = args.Duplicatedata
if single:
    #读取tpm
    a = {}
    a_tpm_file = open(args.A)
    for lines in a_tpm_file.readlines()[1:]:
        line = lines.strip().split()
        a[line[0]] = line[1:]

    b = {}
    b_tpm_file = open(args.B)
    for lines in b_tpm_file.readlines()[1:]:
        line = lines.strip().split()
        b[line[0]] = line[1:]
    #构建CLOUD输入文件格式
    head1 = ['eS1t'+str(i) for i in range(1,10)]
    head2 = ['eS2t'+str(i) for i in range(1,10)]
    head = head1 + head2
    print('\t'.join(head),file=outfile)

    for lines in open(args.T).readlines():
        line = lines.strip().split()
        print('\t'.join(b[line[0]]),'\t'.join(a[line[1]]),sep='\t',file=outfile)

    outfile.close()
    a_tpm_file.close()
    b_tpm_file.close()
elif duplicate:
    # 读取tpm
    a = {}
    a_tpm_file = open(args.A)
    for lines in a_tpm_file.readlines()[1:]:
        line = lines.strip().split()
        a[line[0]] = line[1:]

    b = {}
    b_tpm_file = open(args.B)
    for lines in b_tpm_file.readlines()[1:]:
        line = lines.strip().split()
        b[line[0]] = line[1:]
    # 构建CLOUD输入文件格式
    head1 = ['TPC','TPCA']
    head2 = ['eP'+str(i) for i in range(1,10)]
    head3 = ['eC'+str(i) for i in range(1,10)]
    head4 = ['eA'+str(i) for i in range(1,10)]
    head = head1+head2+head3+head4
    print('\t'.join(head),file=outfile)
    for lines in open(args.T).readlines()[1:]:
        line = lines.strip().split()
        print(line[0],line[1],'\t'.join(b[line[2]]),'\t'.join(b[line[3]]),'\t'.join(a[line[4]]),sep='\t',file=outfile)

    outfile.close()
    a_tpm_file.close()
    b_tpm_file.close()