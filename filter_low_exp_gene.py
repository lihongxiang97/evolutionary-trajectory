#!/usr/bin/env python
# -*- coding: utf-8 -*-
#writed by Lhx on 2/15/2022
import sys
if len(sys.argv) != 4 and len(sys.argv) != 6:
    print("Usage: python filter_low_exp_gene.py dup_gene_file dupgene_species_TPM_file dupgene_species_cutoff_value ")
    print("##################format#################")
    print("dup_gene_file format:dup1\tdup2")
    print("TPM file format : gene_id\tsample1_TPM\tsample2_TPM\tsample3_TPM......")
    print("Cutoff_value is that if both dupgenes' all samples' TPMs are lower than it,the dupgene pair will be filtered!")
    print("#########################################")
    print("It will create outfile:dup_gene_file.filter_low_dup_exp!")
    print("If you want to filter ancestor gene exp lower than cutoff,you can use:")
    print("python filter_low_exp_gene.py dup_gene_file dupgene_species_TPM_file dupgene_species_cutoff_value ancestor_species_TPM_file ancestor_species_cutoff_value")
    print("Now dup_gene_file format should be : dup1\tdup2\tancestor")
    print("Than it will create outfile:dup_gene_file.filter_low_exp!")
    sys.exit()

########################read dup_gene_file###############################
try:
    print("Now reading dup_gene_file......\n")
    dup_gene_f = open(sys.argv[1])
    dup_gene = dup_gene_f.readlines()
    dup_gene_f.close()
except FileNotFoundError:
    print("Cannot open "+sys.argv[1]+"!")
    sys.exit()

########################read TPM_file##############################
try:
    print("Now reading dupgene_species_TPM_file......\n")
    tpm_f = open(sys.argv[2])
    tpm = tpm_f.readlines()
    tpm_f.close()
    tpm_dict = {}
    for lines in tpm:
        if not lines.startswith('Gene_id'):
            line = lines.strip().split()
            tpm_dict[line[0]] = '\t'.join(line[1:])
except FileNotFoundError:
    print("Cannot open "+sys.argv[2]+"!")
    sys.exit()

########################filter dup_gene pairs lower than cutoff#################################
print("Now filtering dup_gene pairs lower than cutoff......\n")
cutoff = float(sys.argv[3])
outfile = open(sys.argv[1]+".filter_low_dup_exp","w")
for lines in dup_gene:
    line = lines.strip().split()
    a = 0
    b = 0
    for t in tpm_dict[line[0]].split('\t'):
        if float(t) >= cutoff:
            a = 1
    for t in tpm_dict[line[1]].split('\t'):
        if float(t) >= cutoff:
            b = 1
    if a == 1 and b == 1:
        print(lines.strip(),file=outfile)
outfile.close()
print("Done!"+sys.argv[1]+".filter_low_dup_exp has been created!\n")
########################filter ancestor gene lower than cutoff################################
if len(sys.argv) ==6:
    try:
        a_tpm_f = open(sys.argv[4])
        a_tpm = a_tpm_f.readlines()
        a_tpm_f.close()
        a_tpm_dict = {}
        a_cutoff = float(sys.argv[5])
        outfile = open(sys.argv[1]+".filter_low_exp","w")
        for lines in a_tpm:
            if not lines.startswith('Gene_id'):
                line = lines.strip().split()
                a_tpm_dict[line[0]] = '\t'.join(line[1:])
        for lines in dup_gene:
            line =lines.strip().split()
            a = 0
            b = 0
            c = 0
            for t in tpm_dict[line[0]].split('\t'):
                if float(t) >= cutoff:
                    a = 1
            for t in tpm_dict[line[1]].split('\t'):
                if float(t) >= cutoff:
                    b = 1
            for t in a_tpm_dict[line[2]].split('\t'):
                if float(t) >= a_cutoff:
                    c = 1
            if a == 1 and b == 1 and c == 1:
                print(lines.strip(),file=outfile)
        print("Done!"+sys.argv[1]+".filter_low_exp has been created!\n")
    except FileNotFoundError:
        print("Cannot open "+sys.argv[4]+"!")
        sys.exit()