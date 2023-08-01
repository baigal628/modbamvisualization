import sys
import pysam
import matplotlib.pyplot as plt
import matplotlib.patches as mplpatches
import math
import numpy
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
import seaborn as sns
import pandas as pd
import os
import argparse
from scipy.signal import savgol_filter

parser = argparse.ArgumentParser(description='Visualizing modified bam files and nucleosome bed files',
                                 usage='python[3+] modbamvisualization.py -b file.bam OR -m file-moddata.txt -r "chrN:startpos-stoppos" [options]')

parser.add_argument('-b', '--bamfile', action='store', dest='b',
                    help='modified bam file, should be indexed. Can be filtered to locus or not. Optimal if name does not include dashes.')
parser.add_argument('-m', '--moddatafile', action='store', dest='m',
                    help='data file from earlier run of this program, contains readnames and modification code for each position')
parser.add_argument('-n', '--nucbed', action='store', dest='n',
                    help='nucleosomes.bed file from cawlr sma')
parser.add_argument('-k', '--pca', action='store_true', dest='k',
                    help='whether to plot PCA and run kmeans clustering (default clusters: 3, specify otherwise with -c)')
parser.add_argument('-p', '--plot', action='store_true', dest='p',
                    help='whether to plot reads with modified positions, will cluster, outputs pdf')
parser.add_argument('-o', '--overlay', action='store_true', dest='o',
                    help='whether to plot reads with modifications overlayed with nucleosome positions, required -n')
parser.add_argument('-x', '--remove', action='store_true', dest='x',
                    help='whether to remove highly/chaotically modified reads. This will only work if you specify a threshold.')
parser.add_argument('-t', '--threshold', action='store', dest='t',
                    help='threshold between 0 and 1 to binarize modifications, mod score above this are called as true. Reccommend looking at score dist of positive and negative controls to determine this. If this is not set, will not binarize for plotting.')
parser.add_argument('-c', '--clusters', action='store', dest='c', default='3',
                    help='number of clusters for plotting, deault 3, reccommend looking at PCA to help determine')
parser.add_argument('-r', '--region', action='store', dest='r',
                    help='region to calculate modification for and plot "chrN:startpos-stoppos"')
args = parser.parse_args()
#k make PCA/cluster
#p plot just mod positions
#o plot mod positions with nucleosomes
#x remove highly/chaotically modified reads
#t threshold (0-255) if specified, this will binarize your modification data for plotting
#r region: "chrN:startpos-stoppos"
#c number of clusters

# print(args)

###going to need to remove highly modified reads using the clusters,

compbase = {'A':'T', 'T':'A', 'C':'G', 'G':'C', 'N':'N'}
def getcomp(seq):
    newseq = []
    for base in seq: newseq.append(compbase[base])
    return ''.join(newseq)#newseq[::-1]

coltocode = {(1.,1.,1.):0, (.8,.8,.8):0, (0., 0., 1.):0, (1., 0., 0.):3, (0.5, 0.5, 1.):0,(0.4, 0.8, 1.):0,(0.4, 0.6, 1.):0, (1., 0.5, 0.5):3}
codetocol = {0:(1.,1.,1.),1:(0.4, 0.8, 1.),2:(0.4,0.6, 1.),3:(1., 0., 0.)}
codeToSimpCode = {0:0, 1:0, 2:0, 3:1}
numToc = {-1.:'grey', 0.:'blue', 1.:'orange', 2.:'green', 3.:'purple', 4.:'red', 5.:'yellow'}


if not args.r:
    raise Exception('must specify region')
if not (args.b or args.m):
    raise Exception('must provide .bam file or moddata.txt file')
if args.o and not args.n:
    raise Exception('if you want to overlay nucleosomes, you must provide a nucleosomes.bed file')

numclust = int(args.c)
chr, temp = args.r.split(':')
qstart, qend = temp.split('-')
qstart, qend = int(qstart), int(qend)
plotrange = qend - qstart
threshold = int(float(args.t)*255) if args.t else -1

readnames, readdir, modmatrix, readstdev = [], [], [], []
fileprefix = ''
if args.m:
    if not os.path.isfile(args.m):
        raise Exception('moddata file does not exist')
    else:
        fileprefix = args.m.split('-')[0]
        for line in open(args.m):
            line = line.split('\t')
            readname, dir = line[0], line[1]
            thismat = [int(x) for x in line[2].rstrip().split(',')]
            knownpos = [x for x in thismat if x != -1]

            readstdev.append((numpy.std(knownpos), sum(knownpos)/len(knownpos)))
            readnames.append(readname)
            readdir.append(dir)
            modmatrix.append(thismat)
        print('done loading modified positions')
elif args.b:
    if not os.path.isfile(args.b):
        raise Exception('bam file does not exist')
    else:
        fileprefix = '.'.join(args.b.split('.')[:-1])
        out = open('-'.join([fileprefix, chr, str(qstart), str(qend)]) + '-moddata.txt', 'w')
        samfile = pysam.AlignmentFile(args.b, "rb")
        for s in samfile.fetch(chr, qstart, qend):
            if not s.is_secondary:
                alignstart, alignend = s.reference_start, s.reference_end
                if alignstart < qstart and alignend > qend:
                    readname = s.query_name
                    print(readname)
                    cigar = s.cigartuples

                    ##modified bases are done on raw sequence, need to then get what they look like on genome using CIGAR string
                    # print(s.modified_bases.keys())
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

                    if s.has_tag('MM'):
                        skippedBase = -1 if s.get_tag('MM').split(',', 2)[0][-1] == '?' else 0
                    elif s.has_tag('Mm'):
                        skippedBase = -1 if s.get_tag('Mm').split(',', 2)[0][-1] == '?' else 0
                    else:
                        pass

                    seq = s.query_sequence
                    seqlen = len(seq)
                    if s.is_reverse:  ###need to get compliment of sequence, but not reverse!!
                        seq = getcomp(seq)
                    seqApos = []
                    c = 0
                    for b in seq:
                        if b == 'A':
                            seqApos.append(c)
                        c += 1

                    c = 0
                    for i in seqApos:
                        if c not in ml:
                            ml.insert(c, (i, skippedBase))
                        elif i != ml[c][0]:
                            ml.insert(c, (i, skippedBase))
                        c += 1

                    ml = dict(ml)

                    ref, quer = 0, 0
                    posOnGenome = []
                    for block in cigar:
                        if block[0] in {0, 7, 8}:  # match, consumes both
                            for i in range(block[1]):
                                if quer in ml: posOnGenome.append([ref + alignstart, ml[quer]])
                                ref += 1
                                quer += 1
                        elif block[0] in {1, 4}:  # consumes query
                            quer += block[1]
                        elif block[0] in {2, 3}:  # consumes reference
                            ref += block[1]
                    readnames.append(readname)
                    dirtowrite = '-' if s.is_reverse else '+'
                    readdir.append(dirtowrite)
                    posdict = dict(posOnGenome)
                    thismat = [-1 for x in range(plotrange)]
                    for i in range(qstart, qend):
                        if i in posdict:
                            thismat[i-qstart] = posdict[i]
                    modmatrix.append(thismat)
                    out.write(readname + '\t' + dirtowrite + '\t' + ','.join([str(x) for x in thismat]) + '\n')
        samfile.close()
        out.close()
        print('done parsing region in bam file')

clusterlabels = []
if args.k or args.p or args.o: ##getting clustering done
    modmatrix = numpy.array(modmatrix)
    pca = PCA(n_components=2)
    pca_features = pca.fit_transform(modmatrix)
    pca_df = pd.DataFrame(data=pca_features, columns=['PC1', 'PC2'])

    kmeans = KMeans(init="k-means++", n_clusters=numclust, n_init=4)
    kmeans.fit(pca_features)
    clusterlabels = kmeans.labels_.astype(float)
    # print(min(clusterlabels), max(clusterlabels))
    # subclust = [[] for i in range(numclust)]
    # for i in range(len(modmatrix)):
    #     subclust[int(clusterlabels[i])].append(modmatrix[i])
    # kmeans = KMeans(init="k-means++", n_clusters=numclust + 1, n_init=4)
    # fig, axs = plt.subplots(1, numclust)
    # for i in range(numclust):
    #     thismat = numpy.array(subclust[i])
    #     pca_features = pca.fit_transform(thismat)
    #     pca_df = pd.DataFrame(data=pca_features, columns=['PC1', 'PC2'])
    #     kmeans.fit(pca_features)
    #     theselabels = kmeans.labels_.astype(float)
    #     axs[i].scatter(pca_features[:, 0], pca_features[:, 1], c=[numToc[x] for x in theselabels], alpha=0.3)
    # plt.savefig('-'.join([fileprefix, chr, str(qstart), str(qend)]) + '-kmeans-cluster-subcluster.png', dpi=600)
    if args.k:
        plt.figure()
        plt.scatter(pca_features[:, 0], pca_features[:, 1], c=[numToc[x] for x in kmeans.labels_.astype(float)], alpha=0.3)
        plt.savefig('-'.join([fileprefix, chr, str(qstart), str(qend)]) + '-kmeans-cluster.png', dpi=600)
    print('done clustering')

readToNucPos = {}
if args.n and args.o:
    if not os.path.isfile(args.n):
        raise Exception('nucleosome positions bed file does not exist')
    else:
        for line in open(args.n):
            if line[:5] != 'track':
                line = line.rstrip().split('\t')
                start = int(line[1])
                readname = line[3]
                exonlen = [int(x) for x in line[10].split(',')]
                exonstarts = [int(x) for x in line[11].split(',')]
                readToNucPos[readname] = []
                for i in range(len(exonstarts)):
                    estart = start + exonstarts[i]
                    eend = estart + exonlen[i]
                    if eend < qstart or estart > qend: continue
                    elif estart < qstart: estart = qstart
                    elif eend > qend: eend = qend
                    readToNucPos[readname].append((estart, eend))
        print('done loading nucleosome positions')

readcolbybp = []
clusteravg = {}
if args.p or args.o:
    for z in range(len(modmatrix)):
        thesemods = modmatrix[z]
        currColor = 0
        lasta, lastpos = 0, 0
        colbybp = [0 for x in range(plotrange)]
        for i in range(qstart, qend):
            modprob = thesemods[i - qstart]
            if i - lasta > 5:
                thesebounds = [i - 3, i + 3]
            elif i - lasta > 3:
                thesebounds = [int(i - ((i - lasta) / 2)), i + 3]
            else:
                thesebounds = [(i - (i - lasta)) + 1, i + 3]

            if modprob == -1:
                currColor = 0  # white
            else:
                lasta = i
                if args.t:
                    if modprob >= threshold:
                        currColor = 3  # red
                    else:
                        if readdir[z] == '-':  # read is reverse
                            currColor = 1
                        else:
                            currColor = 2
                else:
                    currColor = modprob
            for j in range(thesebounds[0], thesebounds[1]):
                if j < qend:
                    colbybp[j - qstart] = currColor  # codetocol[currColor]
        readcolbybp.append(colbybp)

    for i in range(len(readcolbybp)):
        thislabel = clusterlabels[i]
        if thislabel not in clusteravg: clusteravg[thislabel] = [0 for x in range(plotrange)]
        for j in range(plotrange):
            if args.t:
                clusteravg[thislabel][j] += codeToSimpCode[readcolbybp[i][j]]
            else:
                clusteravg[thislabel][j] += readcolbybp[i][j]
        readcolbybp[i] = (thislabel, readcolbybp[i], readnames[i], readdir[i])#, [str(round(x, 3)) for x in readstdev[i]])
    for label in clusteravg:
        if not args.t:
            clusteravg[label] = savgol_filter(clusteravg[label], 15, 1)
        thismax = float(max(clusteravg[label]))
        clusteravg[label] = [x/thismax for x in clusteravg[label]]
    readcolbybp.sort()

    # plt.figure()
    # c = 0
    # totreads = numclust * 4 + 1
    # xaxis = [x+qstart for x in range(plotrange)]
    # ticklocs, ticklabels = [], []
    # for i in range(numclust):
    #     # plt.plot(xaxis, [(totreads - c) + x + 0.2 for x in clusteravg[float(i)]], color='black')
    #     # ticklocs.append(0.5 + totreads - c)
    #     # ticklabels.append('Cluster ' + str(i) + '-nosmoothing')
    #     # c += 1
    #     for j in [15,20]:
    #         for k in [1,2,3]:
    #             smoothed = savgol_filter(clusteravg[float(i)], j, k)
    #             plt.plot(xaxis, [(totreads - c) + x + 0.2 for x in smoothed], color='black')
    #             ticklocs.append(0.5 + totreads - c)
    #             ticklabels.append('Cluster ' + str(i) + '-' + str(j) + '-' + str(k))
    #             c += 1
    # plt.yticks(ticklocs, ticklabels, size='small')
    # plt.savefig('-'.join([fileprefix, chr, str(qstart), str(qend)]) + '-diff-smoothing-on-clusters.png', dpi=600)
    print('done smoothing and calculating cluster averages')

if args.p: #plot mod pos
    figheight = 100 if int(len(readcolbybp)/5) + numclust > 98 else 2+int(len(readcolbybp)/5)+numclust
    plt.figure(figsize = (12, figheight))
    ax = plt.axes(frameon=False)
    box = ax.get_position()
    ax.set_position([box.x0 + (box.width*0.22), box.y0, box.width * 0.85, box.height])
    ticklocs, ticklabels = [], []
    c = 0
    xaxis = [x+qstart for x in range(plotrange)]
    totreads = len(readcolbybp) + (numclust*2) #num clusters
    lastcluster = -2
    for colbybp in readcolbybp:
        if colbybp[0] != lastcluster:
            c += 1
            plt.plot(xaxis, [(totreads-c)+(x*1.6) +0.2 for x in clusteravg[colbybp[0]]], color='black')
            ticklocs.append(0.5 + totreads - c)
            ticklabels.append('Cluster ' + str(int(colbybp[0])))
            c += 1
            lastcluster = colbybp[0]
        mydircode = 2 if colbybp[3] == '+' else 1
        readname = colbybp[2] + ' ' + colbybp[3] #+ ' ' + ' '.join(colbybp[4]) #str(round(colbybp[4], 3))
        colbybp = colbybp[1]
        rectshapes = []
        laststart, lastcol = 0, 0
        for i in range(plotrange):
            if colbybp[i] != lastcol:
                rectshapes.append([laststart+qstart, i-laststart, lastcol, totreads-c])
                laststart, lastcol = i, colbybp[i]
        ticklocs.append(0.5 + totreads - c)
        ticklabels.append(readname)
        c += 1

        lastcol = 0
        for i in range(len(rectshapes)):
            rect = rectshapes[i]
            if args.t:
                if rect[1] <= 5 and (rect[2] == 0 or rect[2] == 1): rect[2] = mydircode
                if rect[2] != 0:
                    myalpha = 1 if rect[2] == 3 else 0.3
                    rectangle1 = mplpatches.Rectangle((rect[0], rect[3]), rect[1], 0.9, facecolor=codetocol[rect[2]],
                                                      linewidth=0, alpha=myalpha)
                    ax.add_patch(rectangle1)
            else:
                if rect[1] <= 5 and (rect[2] == 0 or rect[2] == -1) and i < len(rectshapes)-1:
                    rect[2] = int((rectshapes[i-1][2]+rectshapes[i+1][2])/2)
                myalpha = (rect[2]/365.)+0.3
                rectangle1 = mplpatches.Rectangle((rect[0], rect[3]), rect[1], 0.9, facecolor=(rect[2]/255., 0, 1-(rect[2]/255.)),
                                                  linewidth=0, alpha=myalpha)
                ax.add_patch(rectangle1)

    ax.set_xlim(qstart, qend)
    ax.set_ylim(0,totreads+2)
    plt.yticks(ticklocs, ticklabels, size='small')
    plotname = [fileprefix, chr, str(qstart), str(qend), args.c + 'clusters']
    if args.t:
        plotname.append(str(args.t) + 'threshold')
    plt.savefig('-'.join(plotname) + '-modifiedReads.pdf', format='pdf', dpi=600)
    print('done plotting modified reads')





if args.o: #plot mod pos with overlaid nucleosomes
    figheight = 100 if int(len(readcolbybp)/5) + numclust > 98 else 2+int(len(readcolbybp)/5)+numclust
    plt.figure(figsize = (12, figheight))
    ax = plt.axes(frameon=False)
    box = ax.get_position()
    ax.set_position([box.x0 + (box.width * 0.22), box.y0, box.width * 0.85, box.height])
    ticklocs, ticklabels = [], []
    c = 0
    xaxis = [x+qstart for x in range(plotrange)]
    totreads = len(readcolbybp) + (numclust*2) #num clusters
    lastcluster = -2
    for colbybp in readcolbybp:
        if colbybp[0] != lastcluster:
            c += 1
            plt.plot(xaxis, [(totreads-c)+(x*1.6)+0.2 for x in clusteravg[colbybp[0]]], color='black')
            ticklocs.append(0.5+ totreads-c)
            ticklabels.append('Cluster ' + str(int(colbybp[0])))
            c += 1
            lastcluster = colbybp[0]
        mydircode = 2 if colbybp[3] == '+' else 1
        readname = colbybp[2] + ' ' + colbybp[3]
        colbybp = colbybp[1]
        rectshapes = []
        laststart, lastcol = 0, 0
        for i in range(plotrange):
            if colbybp[i] != lastcol:
                rectshapes.append([laststart+qstart, i-laststart, lastcol, totreads-c])
                laststart, lastcol = i, colbybp[i]

        for block in readToNucPos[readname.split(' ')[0]]:
            rectangle1 = mplpatches.Rectangle((block[0], totreads-c), block[1]-block[0], 0.9, facecolor='grey',
                                              linewidth=0, alpha=0.5)
            ax.add_patch(rectangle1)
        ticklocs.append(0.5 + totreads-c)
        ticklabels.append(readname)
        c += 1

        for i in range(len(rectshapes)):
            rect = rectshapes[i]
            if args.t:
                if rect[1] <= 5 and (rect[2] == 0 or rect[2] == 1): rect[2] = mydircode
                if rect[2] != 0:
                    myalpha = 1 if rect[2] == 3 else 0.3
                    rectangle1 = mplpatches.Rectangle((rect[0], rect[3]+0.25), rect[1], 0.4, facecolor=codetocol[rect[2]],
                                                      linewidth=0, alpha=myalpha)
                    ax.add_patch(rectangle1)
            else:
                if rect[1] <= 5 and (rect[2] == 0 or rect[2] == -1) and i < len(rectshapes)-1:
                    rect[2] = int((rectshapes[i-1][2]+rectshapes[i+1][2])/2)
                myalpha = (rect[2]/365.)+0.3
                rectangle1 = mplpatches.Rectangle((rect[0], rect[3]+0.25), rect[1], 0.4, facecolor=(rect[2]/255., 0, 1-(rect[2]/255.)),
                                                  linewidth=0, alpha=myalpha)
                ax.add_patch(rectangle1)
    ax.set_xlim(qstart, qend)
    ax.set_ylim(0,totreads+2)
    plt.yticks(ticklocs, ticklabels, size='small')
    plotname = [fileprefix, chr, str(qstart), str(qend), args.c + 'clusters']
    if args.t:
        plotname.append(str(args.t) + 'threshold')
    plt.savefig('-'.join(plotname) + '-modifiedReadsWithNucleosomePositions.pdf', format='pdf', dpi=600)
    print('done plotting modified reads overlaid with nucleosomes')
