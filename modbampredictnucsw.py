import sys
import pysam
import matplotlib.pyplot as plt
import matplotlib.patches as mplpatches
import math
import numpy as np
# from sklearn.cluster import KMeans
# from sklearn.decomposition import PCA
# import seaborn as sns
import pandas as pd
import os
import argparse
# from scipy.signal import savgol_filter
# from hmmlearn import hmm
# np.random.seed(42)







def getThreshFromFile(file):
    threshold, lowthresh, posprob = None, None, None
    for line in open(file):
        line = line.rstrip().split('\t')
        if line[0] == 'predicted threshold:': 
            threshold = int(float(line[1])*255)
            t = float(line[1])
        elif line[0] == 'lowest possible threshold:': lowthresh = int(float(line[1])*255)
        elif line[0] == 'pos nuc prob above threshold:': posprob = float(line[1])
    return t, threshold, lowthresh, posprob




def predictOpenChromatin(modmatrix, plotrange, windowsize, threshold, lowthresh, posprob, isbothstrands):
    predictednuc = []
    minacount = int(round(windowsize / 4.)) if not isbothstrands else int(round(windowsize/2.))
    for z in range(len(modmatrix)):
        thesemods = modmatrix[z]
        nucpred = [0 for x in range(plotrange)]
        for i in range(plotrange - windowsize):
            thisspan = thesemods[i:i + windowsize]
            ascores = [x for x in thisspan if x != -1]
            # if len(ascores) > 0: thisscore = sum(ascores)/len(ascores)
            # else: thisscore = 0
            highas = [x for x in ascores if x >= lowthresh]
            # if thisscore >= threshold and len(ascores) >= 6:
            if len(ascores) >= minacount and len(highas) / len(ascores) >= posprob and sum(highas) / len(
                    highas) >= threshold:
                for j in range(windowsize):
                    nucpred[i + j] = 1

        # colbybp = [1 if x >= threshold else 0 for x in thesemods]
        # colbybp = [0 for x in range(plotrange)]
        # for i in range(plotrange):
        #     if thesemods[i] >= threshold: colbybp[i] = 1
        #     elif thesemods[i] >= 0: colbybp[i] = 2
        # readcolbybp.append(colbybp)
        predictednuc.append(nucpred)
    return predictednuc


# threshold = None

# figheight = 100 if int(len(readcolbybp)/5) > 98 else 2+int(len(readcolbybp)/5)
# plt.figure(figsize = (12, figheight))
# ax = plt.axes(frameon=False)
# box = ax.get_position()
# ax.set_position([box.x0 + (box.width*0.22), box.y0, box.width * 0.85, box.height])
# ticklocs, ticklabels = [], []
# xaxis = [x+qstart for x in range(plotrange)]
###num clusters
# lastcluster = -2



def predictNucleosomePos(predictednuc, plotrange, separatenucs):
    c = 0
    # totreads = len(predictednuc)
    newnucpos = []
    for j in range(len(predictednuc)):
        nucpred = predictednuc[j]
        nucshapes = []
        nucstart, nuccol = 0,0
        lastred = 0
        lastblocks = []
        for i in range(plotrange):
            if nucpred[i] != nuccol or i == plotrange-1:
                if nuccol == 0:
                    if separatenucs:
                        blocklen = i-nucstart
                        blockcenter = nucstart + int(blocklen/2)
                        if blocklen/147 >= 0.9 and blocklen/147 <= 2:
                            if len(lastblocks) > 0:
                                for b in lastblocks:
                                    nucshapes.append((b, 147))
                                lastblocks = []
                            blockstart = blockcenter-73
                            lastblocks.append(blockstart)
                        elif blocklen/147 < 0.9 and nucstart-lastred < 5 and len(lastblocks) > 0:
                            nucstart = lastblocks[0]
                            blocklen = i - nucstart
                            blockcenter = nucstart + int(blocklen / 2)
                            lastblocks = []
                        if blocklen/147 > 2:
                            maxnuc = int((blocklen+5)/(147+5))
                            minnuc = int((blocklen+20)/(147+20))
                            bestgap, bestnuc, minerror = 0, 0, 1
                            for l in range(minnuc, maxnuc+1):
                                for m in range(5,21):
                                    error = blocklen/((147*l) + (m * (l-1)))
                                    if error-int(error) < minerror:
                                        bestgap, bestnuc, minerror = m, l, error-int(error)
                            totblocksize = (bestnuc*147) + (bestgap*(bestnuc-1))
                            blockstart = blockcenter - int(totblocksize / 2)
                            if len(lastblocks) > 0:
                                for b in lastblocks:
                                    nucshapes.append((b, 147))
                                lastblocks = []
                            for l in range(bestnuc):
                                lastblocks.append(blockstart)
                                blockstart += 147 + bestgap
                            # for l in range(bestnuc):
                            #     nucshapes.append([blockstart + qstart, 147, 'grey', totreads - c])
                            #     blockstart += 147 + bestgap
                    else:
                        nucshapes.append((nucstart, i-nucstart))
                if nucpred[i] == 1: lastred = i
                nucstart, nuccol = i, nucpred[i]
        if len(lastblocks) > 0:
            for b in lastblocks:
                nucshapes.append((b, 147))
        newnucpos.append(nucshapes)
    return newnucpos


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Predicting nucleosome positions using a sliding window',
                                     usage='python[3+] modbamvisualization.py -m file-moddata.txt -r "chrN:startpos-stoppos" -t predictedThreshold.tsv [options]')

    parser.add_argument('-m', '--moddatafile', action='store', dest='m',
                        help='data file from modbedvisualization, contains readnames and modification code for each position. The region this was run with should be the same as the inputted region.')
    parser.add_argument('-t', '--thresholdfile', action='store', dest='t',
                        help='file generated from predictthreshold.py containing predicted threshold, low threshold, and pos prob. You can make this file yourself as long as you follow the same formatting.')
    parser.add_argument('-r', '--region', action='store', dest='r',
                        help='region to calculate nucleosomes for. make sure this matches the region modbedvisualization was run with. "chrN:startpos-stoppos"')
    parser.add_argument('-w', '--windowsize', action='store', dest='w', default='30',
                        help='window size to scan, larger is more stringent for searching for larger regions of open chromatin. Default: 30, unit is base pairs')
    parser.add_argument('-o', '--predopenonly', action='store_true', dest='o',
                        help='whether to only predict open chromatin and not nucleosomes, default: False')

    args = parser.parse_args()



    args.r = args.r.split(':')
    chr = args.r[0]
    qstart, qend = [int(x) for x in args.r[1].split('-')]
    plotrange = qend-qstart

    t, threshold, lowthresh, posprob = getThreshFromFile(args.t)

    windowsize = int(args.w)

    fileprefix = args.m.split('-')[0]

    readnames, readdir, modmatrix, readstdev = [], [], [], []
    aposlist, ascorelist = [], []
    for line in open(args.m):
        line = line.split('\t')
        readname, dir = line[0], line[1]
        thismat = [int(x) for x in line[2].rstrip().split(',')]
        knownpos = [x for x in thismat if x != -1]

        readstdev.append((np.std(knownpos), sum(knownpos) / len(knownpos)))
        readnames.append(readname)
        readdir.append(dir)

        modmatrix.append(thismat)

    predictednuc = predictOpenChromatin(modmatrix, plotrange, windowsize, threshold, lowthresh, posprob, False)

    print('done predicting open chromatin')

    finalnucpred = predictNucleosomePos(predictednuc, plotrange, not args.o)
    print('done predicting nucleosomes')

    out = open(fileprefix + '_' + str(t) + 'threshold-predictednuc.bed', 'w')

    for i in range(len(finalnucpred)):
        readname = readnames[i]
        nucshapes = finalnucpred[i]
        line = [chr, str(qstart+ nucshapes[0][0]), str(qstart+nucshapes[-1][0]+147), readname, '0', readdir[i], str(qstart+nucshapes[0][0]), str(qstart+nucshapes[-1][0]+147), '0', str(len(nucshapes)), ','.join(['147' for x in range(len(nucshapes))]), ','.join([str(x[0]-nucshapes[0][0]) for x in nucshapes])]
        out.write('\t'.join(line) + '\n')
    out.close()


# c = 0
# totreads = len(readcolbybp)
# out = open(fileprefix + '_' + str(threshold) + 'threshold-predictednuc.bed', 'w')
# for j in range(len(readcolbybp)):
#     colbybp = readcolbybp[j]
#     nucpred = predictednuc[j]
#     rectshapes = []
#     laststart, lastcol = 0, 0
#     nucshapes = []
#     nucstart, nuccol = 0,0
#     lastred = 0
#     lastblocks = []
#     for i in range(plotrange):
#         # if colbybp[i] != lastcol or i == plotrange-1:
#         #     rectshapes.append([laststart+qstart, i-laststart, lastcol, totreads-c])
#         #     laststart, lastcol = i, colbybp[i]
#
#         if nucpred[i] != nuccol or i == plotrange-1:
#             if nuccol == 0:
#                 if not args.o:
#                     blocklen = i-nucstart
#                     blockcenter = nucstart + int(blocklen/2)
#                     # print(nucstart-lastred)
#                     if blocklen/147 >= 0.9 and blocklen/147 <= 2:
#                         if len(lastblocks) > 0:
#                             for b in lastblocks:
#                                 nucshapes.append([b, 147, 'grey', totreads - c])
#                             lastblocks = []
#                         blockstart = blockcenter-73
#                         lastblocks.append(blockstart+qstart)
#                     elif blocklen/147 < 0.9 and nucstart-lastred < 5 and len(lastblocks) > 0:
#                         nucstart = lastblocks[0]-qstart
#                         # print(blocklen, i-nucstart, (i-nucstart)/147)
#                         blocklen = i - nucstart
#                         blockcenter = nucstart + int(blocklen / 2)
#                         lastblocks = []
#                     if blocklen/147 > 2:
#                         maxnuc = int((blocklen+5)/(147+5))
#                         minnuc = int((blocklen+20)/(147+20))
#                         bestgap, bestnuc, minerror = 0, 0, 1
#                         for l in range(minnuc, maxnuc+1):
#                             for m in range(5,21):
#                                 error = blocklen/((147*l) + (m * (l-1)))
#                                 if error-int(error) < minerror:
#                                     bestgap, bestnuc, minerror = m, l, error-int(error)
#                         # print(minerror)
#                         totblocksize = (bestnuc*147) + (bestgap*(bestnuc-1))
#                         blockstart = blockcenter - int(totblocksize / 2)
#                         if len(lastblocks) > 0:
#                             for b in lastblocks:
#                                 nucshapes.append([b, 147, 'grey', totreads - c])
#                             lastblocks = []
#                         for l in range(bestnuc):
#                             lastblocks.append(blockstart + qstart)
#                             blockstart += 147 + bestgap
#                         # for l in range(bestnuc):
#                         #     nucshapes.append([blockstart + qstart, 147, 'grey', totreads - c])
#                         #     blockstart += 147 + bestgap
#                 else:
#                     nucshapes.append([nucstart+qstart, i-nucstart, 'grey', totreads-c])
#                 # print(i-nucstart, (i-nucstart)/147)
#             if nucpred[i] == 1: lastred = i
#             nucstart, nuccol = i, nucpred[i]
#     if len(lastblocks) > 0:
#         for b in lastblocks:
#             nucshapes.append([b, 147, 'grey', totreads - c])
#
#     readname = readnames[j]
#     line = [chr, str(nucshapes[0][0]), str(nucshapes[-1][0]+147), readname, '0', readdir[j], str(nucshapes[0][0]), str(nucshapes[-1][0]+147), '0', str(len(nucshapes)), ','.join(['147' for x in range(len(nucshapes))]), ','.join([str(x[0]-nucshapes[0][0]) for x in nucshapes])]
#     out.write('\t'.join(line) + '\n')
# print('done predicting nucleosomes and writing output')












