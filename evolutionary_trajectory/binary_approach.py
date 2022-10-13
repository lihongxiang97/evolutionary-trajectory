#!/usr/bin/env python
# -*- coding: utf-8 -*-
#writed by Lhx on 2/24/2022
import sys
import matplotlib.pyplot as plt
if len(sys.argv)!=6:
    print("Usage: python binary_approach.py trips_file A_TPM_file B_TPM_file A_TPM_cutoff B_TPM_cutoff")
    sys.exit()
#######################read A_TPM_file##########################
a_tpm={}
try:
    a_tpm_file = open(sys.argv[2])
    for lines in a_tpm_file.readlines():
        line = lines.strip().split()
        a_tpm[line[0]] = []
        for i in line[1:]:
            a_tpm[line[0]].append(i)
    a_tpm_file.close()
except  FileNotFoundError:
    print("Cannot open "+sys.argv[2]+"!")
    sys.exit()
#######################read B_TPM_file##########################
b_tpm={}
try:
    b_tpm_file = open(sys.argv[3])
    for lines in b_tpm_file.readlines():
        line = lines.strip().split()
        b_tpm[line[0]] = []
        for i in line[1:]:
            b_tpm[line[0]].append(i)
    b_tpm_file.close()
except FileNotFoundError:
    print("Cannot open "+sys.argv[3]+"!")
    sys.exit()
########################assess divergence (D) by comparing the tissues in which genes are expressed####################
try:
    trip_file = open(sys.argv[1])
    outfile = open(sys.argv[1]+'.binary_results','w')
    print('Parent\tChild\tAncestor\tD_P,A\tD_C,A\tD_P+C,A\tClassification',file=outfile)
    for lines in trip_file.readlines():
        if not lines.startswith('Parent'):
            line = lines.strip().split()
            D_PA = 0
            for i in range(0,len(a_tpm[line[0]])):
                if (float(a_tpm[line[0]][i]) > float(sys.argv[4]) and float(b_tpm[line[2]][i]) > float(sys.argv[5])) or (float(a_tpm[line[0]][i]) < float(sys.argv[4]) and float(b_tpm[line[2]][i]) < float(sys.argv[5])):
                    pass
                else:
                    D_PA = 1
            D_CA = 0
            for i in range(0,len(a_tpm[line[1]])):
                if (float(a_tpm[line[1]][i]) > float(sys.argv[4]) and float(b_tpm[line[2]][i]) > float(sys.argv[5])) or (float(a_tpm[line[1]][i]) < float(sys.argv[4]) and float(b_tpm[line[2]][i]) < float(sys.argv[5])):
                    pass
                else:
                    D_CA = 1
            D_PCA = 0
            for i in range(0,len(a_tpm[line[0]])):
                if ((float(a_tpm[line[0]][i]) > float(sys.argv[4]) or float(a_tpm[line[1]][i]) > float(sys.argv[4])) and float(b_tpm[line[2]][i]) > float(sys.argv[5])) or ((float(a_tpm[line[0]][i]) < float(sys.argv[4]) and float(a_tpm[line[1]][i]) < float(sys.argv[4])) and float(b_tpm[line[2]][i]) < float(sys.argv[5])):
                    pass
                else:
                    D_PCA = 1
            if D_PA == 0 and D_CA == 0:
                print(lines.strip(),D_PA,D_CA,D_PCA,'Conservation',sep='\t',file=outfile)
            if D_PA == 1 and D_CA == 0:
                print(lines.strip(),D_PA,D_CA,D_PCA,'Neofunctionalization(Parent)',sep='\t',file=outfile)
            if D_PA == 0 and D_CA == 1:
                print(lines.strip(),D_PA,D_CA,D_PCA,'Neofunctionalization(Child)',sep='\t',file=outfile)
            if D_PA == 1 and D_CA == 1 and D_PCA == 0:
                print(lines.strip(),D_PA,D_CA,D_PCA,'Subfunctionalization',sep='\t',file=outfile)
            if D_PA == 1 and D_CA == 1 and D_PCA == 1:
                print(lines.strip(),D_PA,D_CA,D_PCA,'Specialization',sep='\t',file=outfile)
    trip_file.close()
    outfile.close()

except FileNotFoundError:
    print("Cannot open "+sys.argv[1]+"!")
    sys.exit()
#############################count#####################################
binary_results = open(sys.argv[1]+'.binary_results')
outfile = open(sys.argv[1]+'.binary_results_number','w')
lines = binary_results.read()
con = lines.count('Conservation')
pneo = lines.count('Neofunctionalization(Parent)')
cneo = lines.count('Neofunctionalization(Child)')
sub = lines.count('Subfunctionalization')
spe = lines.count('Specialization')
print('Classification\tnumber',file=outfile)
print('Conservation\t',con,file=outfile)
print('Neofunctionalization(Parent)\t',pneo,file=outfile)
print('Neofunctionalization(Child)\t',cneo,file=outfile)
print('Subfunctionalization\t',sub,file=outfile)
print('Specialization\t',spe,file=outfile)
###############draw#######################
data = [con,pneo,cneo,sub,spe]
labels = ['Conservation','Neofunctionalization(Parent)','Neofunctionalization(Child)','Subfunctionalization','Specialization']
fig = plt.figure(0)
plt.bar(range(len(data)),data,tick_label=labels)
plt.xticks(rotation=-15)
plt.tight_layout()
plt.savefig(sys.argv[1]+'.binary_results_number.png',format='png')
plt.close(0)
binary_results.close()
outfile.close()
