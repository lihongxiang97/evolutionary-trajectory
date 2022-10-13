#!/usr/bin/env python
# -*- coding: utf-8 -*-
#writed by Lhx on 4/5/2022
import sys
if len(sys.argv) != 3:
    print('Usage:python script.py classification bed12')
else:
    d={}
    with open(sys.argv[2]) as f:
        for lines in f:
            line = lines.strip().split()
            d[line[3]] = line[9]
    l=[]
    out = open('DNA_or_RNA_mediated.classifications','w')
    with open(sys.argv[1]) as f:
        for lines in f:
            line = lines.strip().split()
            if int(d[line[0]]) > 1 and int(d[line[1]]) > 1:
                print(lines.strip(), 'DNA-mediated', sep='\t', file=out)
                l.append((line[3], 'DNA-mediated'))
            elif int(d[line[0]]) > 1 and int(d[line[1]]) == 1:
                print(lines.strip(), 'RNA-mediated', sep='\t', file=out)
                l.append((line[3], 'RNA-mediated'))
            elif int(d[line[0]]) == 1 and int(d[line[1]]) == 1:
                print(lines.strip(), 'Unknown', sep='\t', file=out)
                l.append((line[3], 'Unknown'))
    out.close()
    count = open('DNA_or_RNA_mediated.count','w')
    print('Conservation','DNA-mediated',l.count(('Conservation','DNA-mediated')),sep='\t',file=count)
    print('Conservation', 'RNA-mediated', l.count(('Conservation', 'RNA-mediated')), sep='\t', file=count)
    print('Neofunctionalization(Parent)', 'DNA-mediated', l.count(('Neofunctionalization(Parent)', 'DNA-mediated')), sep='\t', file=count)
    print('Neofunctionalization(Parent)', 'RNA-mediated', l.count(('Neofunctionalization(Parent)', 'RNA-mediated')), sep='\t', file=count)
    print('Neofunctionalization(Child)', 'DNA-mediated', l.count(('Neofunctionalization(Child)', 'DNA-mediated')), sep='\t', file=count)
    print('Neofunctionalization(Child)', 'RNA-mediated', l.count(('Neofunctionalization(Child)', 'RNA-mediated')), sep='\t', file=count)
    print('Specialization', 'DNA-mediated', l.count(('Specialization', 'DNA-mediated')), sep='\t', file=count)
    print('Specialization', 'RNA-mediated', l.count(('Specialization', 'RNA-mediated')), sep='\t', file=count)
    print('Subfunctionalization', 'DNA-mediated', l.count(('Subfunctionalization', 'DNA-mediated')), sep='\t', file=count)
    print('Subfunctionalization', 'RNA-mediated', l.count(('Subfunctionalization', 'RNA-mediated')), sep='\t', file=count)
    count.close()