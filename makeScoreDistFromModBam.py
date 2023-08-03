import sys
import pysam
import matplotlib.pyplot as plt
import matplotlib.patches as mplpatches
import math
import numpy as np



###usage: python3 makeScoreDistFromModBam.py modtype(ex 6mA) file1.bam file2.bam (file3.bam etc)
##make sure bam files have indexes

###this takes every position in the file that has a modification score and includes it in the final curve

typesOfMods = {'5mC':[('C', 0, 'm')], '5hmC': [('C', 0, 'h')], '5fC': [('C', 0, 'f')], '5caC': [('C', 0, 'c')],
               '5hmU': [('T', 0, 'g')], '5fU': [('T', 0, 'e')], '5caU': [('T', 0, 'b')],
               '6mA': [('A', 0, 'a'), ('A', 0, 'Y')], '8oxoG': [('G', 0, 'o')], 'Xao': [('N', 0, 'n')]}

allmodcount = []
outname = []
for file in range(len(sys.argv)-1):
    outname.append('.'.join(sys.argv[file+1].split('.')[:-1]))
    modCount = [0 for x in range(255)]
    allmods, lowmods, highmods = [], [], []
    samfile = pysam.AlignmentFile(sys.argv[file+1],"rb")
    for s in samfile:
        if not s.is_secondary and s.is_mapped:
            posstag = typesOfMods[args.y]
            if s.is_reverse: posstag = [(x[0], 1, x[2]) for x in posstag]
            ml = None
            for t in posstag:
                if t in s.modified_bases:
                    ml = s.modified_bases[t]
                    break
            if not ml:
                continue
            # print(s.modified_bases.keys())
            # print(ml)
            for i in ml:
                modCount[i[1]-1] += 1
                allmods.append(i[1])
                if i[1] <= 127: lowmods.append(i[1])
                else: highmods.append(i[1])
    allmodcount.append(modCount)
    samfile.close()
    print('done processing ' + sys.argv[file+1])

plt.figure()
filenamesforout = ''
for i in range(len(allmodcount)):
    tot = float(sum(allmodcount[i]))
    allmodcount[i] = [x/tot for x in allmodcount[i]]
    plt.plot([x/255 for x in list(range(255))], allmodcount[i], label = '.'.join(sys.argv[i+1].split('.')[:-1]))
    filenamesforout += '.'.join(sys.argv[i+1].split('.')[:-1]) + '-'
plt.legend()
plt.savefig('-'.join(outname) + '-scoreDist.png', dpi=600)

