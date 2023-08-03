import sys
import pysam
import matplotlib.pyplot as plt
import matplotlib.patches as mplpatches
import math
import numpy as np
import argparse

parser = argparse.ArgumentParser(description='Predicting the best modification threshold and dist from positive and negative control bam files. Will output a tsv file containing those values.',
                                 usage='python[3+] predictthreshold.py -p pos.bam -n neg.bam OR JUST -s sample.bam')
parser.add_argument('-p', '--posbam', action='store', dest='p',
                    help='positive control bam file, should be indexed. Can be filtered to locus or not. Should be nucleosome free and modified the same amount as your nuclear DNA')
parser.add_argument('-n', '--negbam', action='store', dest='n',
                    help='negative control bam file, should be indexed. Can be filtered to locus or not. Should be nucleosome free and unmodified')
parser.add_argument('-s', '--samplebam', action='store', dest='s',
                    help='sample bam file, should be indexed. Can be filtered to locus or not. ONLY USE THIS OPTION IF YOU HAVE NO CONTROLS, SUBOPTIMAL')
# parser.add_argument('-e', '--extra', action='store', dest='e',
#                     help='sample bam file, should be indexed. Can be filtered to locus or not. ONLY USE THIS OPTION IF YOU HAVE NO CONTROLS, SUBOPTIMAL')
parser.add_argument('-o', '--output', action='store', dest='o', default='predictedThreshold',
                    help='output file prefix, default is predictedThreshold')
args = parser.parse_args()


def getmoddist(filename):
    modCount = [0 for x in range(255)]
    samfile = pysam.AlignmentFile(filename, "rb")
    for s in samfile:
        if not s.is_secondary and s.is_mapped:
            if s.is_reverse:
                if ('A', 1, 'a') in s.modified_bases:
                    ml = s.modified_bases[('A', 1, 'a')]
                elif ('A', 1, 'Y') in s.modified_bases:
                    ml = s.modified_bases[('A', 1, 'Y')]
                else: continue
            else:
                if ('A', 0, 'a') in s.modified_bases: ml = s.modified_bases[('A', 0, 'a')]
                if ('A', 0, 'Y') in s.modified_bases:
                    ml = s.modified_bases[('A', 0, 'Y')]
                else: continue
            for i in ml:
                modCount[i[1] - 1] += 1
    samfile.close()
    print('done processing ' + filename)
    tot = float(sum(modCount))
    modCount =  [x/tot for x in modCount]
    return(modCount)

#pos prob > 2x neg prob
#R10: 218
#R9: 183

# >=x% of all pos prob are here
#R10:   50%:83      30%:148
#R9:    50%:218     30%:235

# pos is x fraction of pos+neg prob above i
#R10:   .7:145  .8:218  .9:245
#R9:    .7:131  .8:152  .9:175


if args.p and args.n:
    out = open(args.o + '.tsv', 'w')
    posmodcount = getmoddist(args.p)
    negmodcount = getmoddist(args.n)
    possum, rightsum, possum30, posprob2x = sum(posmodcount), 0, 0, 0
    for i in range(254, -1, -1):
        rightsum += posmodcount[i]
        if rightsum >= possum/2:
            possum30 = i
            break
    for i in range(255):
        if negmodcount[i] > 0 and posmodcount[i]/negmodcount[i] > 2:
            posprob2x = i
            break
    cutoffval = int((possum30+posprob2x)/2)
    lowval = min([possum30, posprob2x]) if min([possum30, posprob2x]) >= 128 else 128
    # print(round(cutoffval/255, 3))
    # print(sum(allmodcount[0][cutoffval:])/sum(allmodcount[0]))
    out.write('predicted threshold:\t' + str(round(cutoffval/255, 3)) + '\n')
    out.write('lowest possible threshold:\t' + str(round(lowval/255, 3)) + '\n')
    out.write('pos nuc prob above threshold:\t' + str(round(sum(posmodcount[cutoffval:])/sum(posmodcount), 3)) + '\n')

elif args.s:
    out = open(args.o + '.tsv', 'w')
    samplemodcount = getmoddist(args.s)
    # sample2 = getmoddist(args.e)
    # fracidiff = {x/100.:[0,0] for x in range(10, 21)}
    samplesum, rightsum = sum(samplemodcount), 0
    # s2sum, rightsum2 = sum(sample2), 0
    cutoffval = 127
    for i in range(254, 127, -1):
        rightsum += samplemodcount[i]
        if rightsum/samplesum >= 0.15:
            cutoffval = i
            break
    out.write('predicted threshold:\t' + str(round(cutoffval / 255, 3)) + '\n')
    out.write('pos nuc prob above threshold:\t' + str(0.5) + '\n')
    #     rightsum2 += sample2[i]
    #     # if rightsum/samplesum < 0.3 and rightsum/samplesum > 0.1: print(i, rightsum/samplesum, rightsum2/s2sum, abs((rightsum/samplesum)-(rightsum2/s2sum)))
    #     rightfrac, r2frac = round(rightsum/samplesum, 2), round(rightsum2/s2sum, 2)
    #     if rightfrac in fracidiff and fracidiff[rightfrac][0] == 0: fracidiff[rightfrac][0] = i
    #     if r2frac in fracidiff and fracidiff[r2frac][1] == 0: fracidiff[r2frac][1] = i
    # for i in range(10, 21):
    #     i = i/100.
    #     print(i, fracidiff[i], abs(fracidiff[i][1]-fracidiff[i][0]))
else:
    print('did not provide enough inputs. Did you provide -n and -p OR -s?')