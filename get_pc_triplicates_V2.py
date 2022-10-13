import sys
from collections import defaultdict
import argparse
parser = argparse.ArgumentParser(description="Identify parent or child gene in duplication and ancestor of the duplication! The input files are from Dupgen_finder!")
parser.add_argument("-A","--a",help="Target species abbrev")
parser.add_argument("-B","--b",help="Outgroup species abbrev")
parser.add_argument("-O","--order",help="If the order of collinearity is A B,input f;else,input r")
args = parser.parse_args()
a_abbrev = args.a
b_abbrev = args.b
order = args.order
#if len(sys.argv)!=4:
#    print("Usage: python scrtpt.py -A A-abbrev(target) -B B-abbrev(outgroup) -O order(f or r)")
#    print("If the target geneid located in col1, please use 'f'")
#    sys.exit()
#############read A_B.collinearity#########################
orp = defaultdict(dict)
try:
    collinearity_file=open(a_abbrev+"_"+b_abbrev+".collinearity",'r')

    for lines in collinearity_file.readlines():
        if not lines.startswith('#'):
            a=lines.strip().split('\t')
            if order=='f':
                orp[a[1]][a[2]]='A'
            elif order=='r':
                orp[a[2]][a[1]]='A'
    collinearity_file.close()

except FileNotFoundError:
    print("Cannot open "+a_abbrev+"_"+b_abbrev+".collinearity!")
    sys.exit()

n=0
orpn={}
orpn1={}
for key1 in orp.keys():
    b=''
    for key2 in orp[key1].keys():
        n+=1
        b=key2
    orpn[key1]=n

    if n==1:
        orpn1[key1]=b
    n=0

##################read A_B.blast file#########################
bla = {}
try:
    blast_file=open(a_abbrev+"_"+b_abbrev+".blast",'r')

    for lines in blast_file.readlines():
        a=lines.strip().split('\t')
        key=a[0] + '\t'+a[1]
        if key not in bla:
            bla[key] = a[11]
        else:
            if float(a[11]) > float(bla[key]):
                bla[key] = a[11] #get the largest blast values
    blast_file.close()

except FileNotFoundError:
    print("Cannot open "+a_abbrev+"_"+b_abbrev+".blast!")
    sys.exit()

blap = {}
temp = {}
for key in bla.keys():
    a=key.split('\t')
    if a[0] not in blap:
        blap[a[0]] = a[1]
        temp[a[0]] = bla[key]
    else:
        if float(bla[key]) > float(temp[a[0]]):
            blap[a[0]] = a[1] #find the best hit gene in outgroup species
            temp[a[0]] = bla[key]

#############searching the 'A A B' triplicates for different modes of duplicated genes#########################           
################WGD######################
try:
    wgd_pairs=open(a_abbrev+".wgd.pairs")
    wgd_out=open(a_abbrev+"_"+b_abbrev+".wgd.trips","w")

    print("Duplicate1","Duplicate2","Ancestral",sep='\t',file=wgd_out)

    for lines in wgd_pairs.readlines():
        if not lines.startswith("Duplicate"):
            a=lines.strip().split('\t')
            if a[0] in orpn and a[2] in orpn:#judge whether dup1 and dup2 syntenic
                if orpn[a[0]] == 1 and orpn[a[2]] == 1:#judge whether the syntenic gene numbers are 1
                    if orpn1[a[0]] == orpn1[a[2]]:#judge whether the only syntenic gene is same
                        print(a[0],a[2],orpn1[a[0]],sep='\t',file=wgd_out)
    wgd_pairs.close()
    wgd_out.close()

except FileNotFoundError:
    print('wgd.pairs not exists!')

################Tandem######################
try:
    tandem_pairs=open(a_abbrev+".tandem.pairs")
    tandem_out=open(a_abbrev+"_"+b_abbrev+".tandem.trips","w")
    tandem_pc_out=open(a_abbrev+"_"+b_abbrev+".tandem_pc.trips","w")

    print("Duplicate1","Duplicate2","Ancestral",sep='\t',file=tandem_out)
    print("Parent","Child","Ancestral",sep='\t',file=tandem_pc_out)

    for lines in tandem_pairs.readlines():
        if not lines.startswith("Duplicate"):
            a=lines.strip().split('\t')
            x=0
            y=0
            if a[0] in blap and a[2] in blap:
                if blap[a[0]] == blap[a[2]]:
                    if a[0] in orp and blap[a[0]] in orp[a[0]]:
                        x=1
                    if a[2] in orp and blap[a[0]] in orp[a[2]]:
                        y=1
            if x==1 or y==1:
                print(a[0], a[2], blap[a[0]], sep='\t', file=tandem_out)
            if x==1 and y==0:
                print(a[0], a[2], blap[a[0]], sep='\t', file=tandem_pc_out)
            if x==0 and y==1:
                print(a[2], a[0], blap[a[0]], sep='\t', file=tandem_pc_out)
#分parent和child gene 原理：parent gene和Ancestral gene有共线性关系，而child gene则没有
    tandem_pairs.close()
    tandem_out.close()
    tandem_pc_out.close()

except FileNotFoundError:
    print('tandem.pairs not exists!')

################Proximal######################
try:
    proximal_pairs=open(a_abbrev+".proximal.pairs")
    proximal_out=open(a_abbrev+"_"+b_abbrev+".proximal.trips","w")
    proximal_pc_out=open(a_abbrev+"_"+b_abbrev+".proximal_pc.trips","w")

    print("Duplicate1", "Duplicate2", "Ancestral", sep='\t', file=proximal_out)
    print("Parent", "Child", "Ancestral", sep='\t', file=proximal_pc_out)

    for lines in proximal_pairs.readlines():
        if not lines.startswith("Duplicate"):
            a = lines.strip().split('\t')
            x=0
            y=0
            if a[0] in blap and a[2] in blap:
                if blap[a[0]] == blap[a[2]]:
                    if a[0] in orp and blap[a[0]] in orp[a[0]]:
                        x=1
                    if a[2] in orp and blap[a[0]] in orp[a[2]]:
                        y=1
            if x==1 or y==1:
                print(a[0], a[2], blap[a[0]], sep='\t', file=proximal_out)
            if x==1 and y==0:
                print(a[0], a[2], blap[a[0]], sep='\t', file=proximal_pc_out)
            if x==0 and y==1:
                print(a[2], a[0], blap[a[0]], sep='\t', file=proximal_pc_out)

    proximal_pairs.close()
    proximal_out.close()
    proximal_pc_out.close()

except FileNotFoundError:
    print('proximal.pairs not exists!')

################Transposed######################
try:
    transposed_pairs=open(a_abbrev+".transposed.pairs")
    transposed_out=open(a_abbrev+"_"+b_abbrev+".transposed.trips","w")
    transposed_pc_out=open(a_abbrev+"_"+b_abbrev+".transposed_pc.trips","w")

    print("Duplicate1", "Duplicate2", "Ancestral", sep='\t', file=transposed_out)
    print("Parent", "Child", "Ancestral", sep='\t', file=transposed_pc_out)

    for lines in transposed_pairs.readlines():
        if not lines.startswith("Duplicate"):
            a = lines.strip().split('\t')
            x = 0
            y = 0
            if a[0] in blap and a[2] in blap:
                if blap[a[0]] == blap[a[2]]:
                    if a[0] in orp and blap[a[0]] in orp[a[0]]:
                        x = 1
                    if a[2] in orp and blap[a[0]] in orp[a[2]]:
                        y = 1
            if x==1 or y==1:
                print(a[0], a[2], blap[a[0]], sep='\t', file=transposed_out)
            if x==1 and y==0:
                print(a[0], a[2], blap[a[0]], sep='\t', file=transposed_pc_out)
            if x==0 and y==1:
                print(a[2], a[0], blap[a[0]], sep='\t', file=transposed_pc_out)

    transposed_pairs.close()
    transposed_out.close()
    transposed_pc_out.close()

except FileNotFoundError:
    print('proximal.pairs not exists!')

################Dispersed######################
try:
    dispersed_pairs = open(a_abbrev + ".dispersed.pairs")
    dispersed_out = open(a_abbrev + "_" + b_abbrev + ".dispersed.trips", "w")
    dispersed_pc_out=open(a_abbrev+"_"+b_abbrev+".dispersed_pc.trips","w")

    print("Duplicate1", "Duplicate2", "Ancestral", sep='\t', file=dispersed_out)
    print("Parent", "Child", "Ancestral", sep='\t', file=dispersed_pc_out)

    for lines in dispersed_pairs.readlines():
        if not lines.startswith("Duplicate"):
            a = lines.strip().split('\t')
            x = 0
            y = 0
            if a[0] in blap and a[2] in blap:
                if blap[a[0]] == blap[a[2]]:
                    if a[0] in orp and blap[a[0]] in orp[a[0]]:
                        x = 1
                    if a[2] in orp and blap[a[0]] in orp[a[2]]:
                        y = 1
            if x==1 or y==1:
                print(a[0], a[2], blap[a[0]], sep='\t', file=dispersed_out)
            if x==1 and y==0:
                print(a[0], a[2], blap[a[0]], sep='\t', file=dispersed_pc_out)
            if x==0 and y==1:
                print(a[2], a[0], blap[a[0]], sep='\t', file=dispersed_pc_out)

    dispersed_pairs.close()
    dispersed_out.close()
    dispersed_pc_out.close()
    
except FileNotFoundError:
    print('proximal.pairs not exists!')

################Singletons######################
singletons={}
try:
    singletons_f1 = open(b_abbrev+".singletons")
    for lines in singletons_f1.readlines():
        if not lines.startswith("GeneID"):
            a=lines.strip().split('\t')
            singletons[a[0]] = 'B'
    singletons_f1.close()

except FileNotFoundError:
    print(b_abbrev+".singletons not exists!")

try:
    singletons_f2 = open(a_abbrev+".singletons")
    singletons_out = open(a_abbrev+"_"+b_abbrev+".singletons",'w')

    print(a_abbrev,b_abbrev,sep='\t',file=singletons_out)

    for lines in singletons_f2.readlines():
        if not lines.startswith("GeneID"):
            a = lines.strip().split('\t')
            if a[0] in orp:
                for key in orp[a[0]].keys():
                    if key in singletons:
                        print(a[0],key,sep='\t',file=singletons_out)
    singletons_f2.close()
    singletons_out.close()

except FileNotFoundError:
    print(a_abbrev+'.singletons not exists!')



    
