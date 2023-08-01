import sys
import pysam
import matplotlib.pyplot as plt
import matplotlib.patches as mplpatches
import math



###usage: python3 makeScoreDistFromModBam.py file1.bam file2.bam (file3.bam etc)
##make sure bam files have indexes

###this takes every position in the file that has a modification score and includes it in the final curve

allmodcount = []
outname = []
for file in range(len(sys.argv)-1):
    outname.append('.'.join(sys.argv[file+1].split('.')[:-1]))
    modCount = [0 for x in range(255)]
    samfile = pysam.AlignmentFile(sys.argv[file+1],"rb")
    for s in samfile:
        if not s.is_secondary and s.is_mapped:
            if s.is_reverse:
                if ('A', 1, 'a') in s.modified_bases:
                    ml = s.modified_bases[('A', 1, 'a')]
                elif ('A', 1, 'Y') in s.modified_bases:
                    ml = s.modified_bases[('A', 1, 'Y')]
                else:
                    pass
            else:
                if ('A', 0, 'a') in s.modified_bases: ml = s.modified_bases[('A', 0, 'a')]
                if ('A', 0, 'Y') in s.modified_bases:
                    ml = s.modified_bases[('A', 0, 'Y')]
                else:
                    pass
            # print(s.modified_bases.keys())
            # print(ml)
            for i in ml:
                modCount[i[1]-1] += 1
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