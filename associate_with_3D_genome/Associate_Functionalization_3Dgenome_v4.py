#! /usr/bin/env python
#Writed by LiHongxiang on 2/28/2022
#修改了compartment的输出，将不存在输出为NA
#edited on 5/2/2022
#修改了输入文件格式，不再是CDROM的输出格式，而是parent\tchild\tancestor\tclassifications
import argparse
parser = argparse.ArgumentParser(description="Associate functionalizalion and 3D genome struction:A/B compartment,TAD,LOOP!")
parser.add_argument("-F",required=True,help="input file, format:parent\tchild\tancestor\tclassifications")
parser.add_argument("-G","--gff",required=True,help="gff,get by 'prepare_gff-v2.pl'")
parser.add_argument("-C","--compartment",help="compartment_file,format:Chr04   1       50000   -0.033  B")
parser.add_argument("-T","--tad",help="TAD,format:Chr10 12195000 12535000")
parser.add_argument("-L","--loop",help="loop,format:Chr01   2500    Chr01   7500")
parser.add_argument("-A","--abbrev",required=True,help="abbrev of generated file")
args = parser.parse_args()
func = args.F
gff = args.gff
com = args.compartment
tad = args.tad
loop = args.loop
abbrev = args.abbrev

##################read input file#########################
d_fcl = {}
func_f = open(func)
print("\033[36mStart reading input file......\n\033[0m")
for lines in func_f.readlines():
    d_tmp = {}
    if not lines.startswith("Parent"):
        line = lines.strip().split()
        d_tmp = {'Parent':line[0],'Child':line[1]}
        if line[3] in d_fcl:
            d_fcl[line[3]].append(d_tmp)
        else:
            d_fcl[line[3]] = [d_tmp]
func_f.close()

######################read gff#############################
d_gff = {}
gff_f = open(gff)
print("\033[36mStart reading gff file......\n\033[0m")
for lines in gff_f.readlines():
    line = lines.strip().split()
    if '-' in line[0]:
        chr = line[0].split('-')[1]
    else:
        chr = line[0]
    d_gff[line[1]] = [chr,line[2],line[3]]
gff_f.close()

######################annotate compartment of duplication######################
if com:
    print("\033[36mStart annotate compartment of duplication......\n\033[0m")
    com_f = open(com)
    com_lines = com_f.readlines()
    com_f.close()
    compartment_c = open(abbrev+"_functionalization_compartment",'w')
    print("Parent\tP_compartment\tChild\tC_compartment\tClassification\tDup_type",file=compartment_c)
    for key,value in d_fcl.items():
        for ids in value:
            x=''
            y=''
            for lines in com_lines:
                line = lines.strip().split()
                if d_gff[ids["Parent"]][0] == line[0]:
                    if int(line[1]) <= int(d_gff[ids["Parent"]][1]) <= int(line[2]):
                        x = line[4]
                if d_gff[ids["Child"]][0] == line[0]:
                    if int(line[1]) <= int(d_gff[ids["Child"]][1]) <= int(line[2]):
                        y = line[4]
            if x and y:
                print(ids["Parent"], x, ids["Child"], y, key, abbrev, sep = '\t', file = compartment_c)
            if x and not y:
                print(ids["Parent"], x, ids["Child"], 'NA', key, abbrev, sep = '\t', file = compartment_c)
            if not x and y:
                print(ids["Parent"], 'NA', ids["Child"], y, key, abbrev, sep = '\t', file = compartment_c)
            if not x and not y:
                print(ids["Parent"], 'NA', ids["Child"], 'NA', key, abbrev, sep = '\t', file = compartment_c)

    compartment_c.close()
    print("\033[34mDone!Generate file:"+abbrev+"_functionalization_compartment!\n\033[0m")
else:
    print("\033[31mIf you want to annotate compartment, please use -C to input compartment file!\n\033[0m")
    pass

#######################annotate TAD of duplication################################
if tad:
    print("\033[36mStart annotate TAD of duplication......\n\033[0m")
    tad_f = open(tad)
    tad_lines = tad_f.readlines()
    tad_f.close()
    tad_c = open(abbrev+"_functionalization_TAD",'w')
    print("Parent\tTAD_name\tChild\tTAD_name\tclassification\tDup_type",file=tad_c)

    for key,value in d_fcl.items():
        for ids in value:
            x=[]
            y=[]

            for lines in tad_lines:
                TAD_name = "TAD"+str(tad_lines.index(lines)+1)
                line = lines.strip().split()
                if d_gff[ids["Parent"]][0] == line[0]:
                #染色体对齐
                    if int(d_gff[ids["Parent"]][1]) <= int(line[2]) and int(line[1]) <= int(d_gff[ids["Parent"]][2]):
                    #TAD区域和gene区域相交
                        x.append(TAD_name)

                if d_gff[ids["Child"]][0] == line[0]:
                    if int(d_gff[ids["Child"]][1]) <= int(line[2]) and int(line[1]) <= int(d_gff[ids["Child"]][2]):
                        y.append(TAD_name)
            if x and y:
                print(ids["Parent"], ','.join(x), ids["Child"], ','.join(y), key, abbrev, sep='\t', file=tad_c)
            if x and not y:
                print(ids["Parent"], ','.join(x), ids["Child"], 'out_TAD', key, abbrev, sep='\t', file=tad_c)
            if not x and y:
                print(ids["Parent"], 'out_TAD', ids["Child"], ','.join(y), key, abbrev, sep='\t', file=tad_c)
            if not x and not y:
                print(ids["Parent"], 'out_TAD', ids["Child"], 'out_TAD', key, abbrev, sep='\t', file=tad_c)
    tad_c.close()
    print("\033[34mDone!Generate file:"+abbrev+"_functionalization_TAD!\n\033[0m")
else:
    print("\033[31mIf you want to annotate TAD, please use -T to input TAD file!\n\033[0m")
    pass
###################annotate loop of duplication###############################
if loop:
    print("\033[36mStart annotate loop of duplication......\n\033[0m")
    loop_f = open(loop)
    loop_lines = loop_f.readlines()
    loop_f.close()
    loop_c = open(abbrev+"_functionalization_loop",'w')
    print("Parent\tLoop_number\tChild\tLoop_number\tclassification\tDup_type",file=loop_c)

    for key,value in d_fcl.items():
        for ids in value:
            x=0
            y=0
            m=[]
            n=[]
            for lines in loop_lines:
                line = lines.strip().split()
                if d_gff[ids["Parent"]][0] == line[0] or d_gff[ids["Parent"]][0] == line[2]:
                    if (int(line[1]) - 2500 <= int(d_gff[ids["Parent"]][2]) and int(line[1]) + 2500 >= int(d_gff[ids["Parent"]][1])) or (int(line[3]) - 2500 <= int(d_gff[ids["Parent"]][2]) and int(line[3]) + 2500 >= int(d_gff[ids["Parent"]][1])):
                        x+=1
                        m.append(line[0:4])
                if d_gff[ids["Child"]][0] == line[0] or d_gff[ids["Child"]][0] == line[2]:
                    if (int(line[1]) - 2500 <= int(d_gff[ids["Child"]][2]) and int(line[1]) + 2500 >= int(d_gff[ids["Child"]][1])) or (int(line[3]) - 2500 <= int(d_gff[ids["Child"]][2]) and int(line[3]) + 2500 >= int(d_gff[ids["Child"]][1])):
                        y+=1
                        n.append(line[0:4])
            if m and n:
                print(ids["Parent"], x, ids["Child"], y, key, abbrev,"Parent_loop:",'\t'.join(j for i in m for j in i),"Child_loop:",'\t'.join(j for i in n for j in i),sep='\t',file=loop_c)
            if m and not n:
                print(ids["Parent"], x, ids["Child"], y, key, abbrev, "Parent_loop:", '\t'.join(j for i in m for j in i),"Child_loop:None", sep='\t', file=loop_c)
            if not m and n:
                print(ids["Parent"], x, ids["Child"], y, key, abbrev, "Parent_loop:None","Child_loop:", '\t'.join(j for i in n for j in i), sep='\t', file=loop_c)
            if not m and not n:
                print(ids["Parent"], x, ids["Child"], y, key, abbrev, "Parent_loop:None","Child_loop:None", sep='\t', file=loop_c)
    loop_c.close()
    print("\033[34mDone!Generate file:"+abbrev+"_functionalization_loop!\n\033[0m")
else:
    print("\033[31mIf you want to annotate loop, please use -L to input loop file!\n\033[0m")
    pass