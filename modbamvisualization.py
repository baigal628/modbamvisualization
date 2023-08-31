import sys
import pysam
import matplotlib.pyplot as plt
import matplotlib.patches as mplpatches
import math
import numpy as np
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
import seaborn as sns
import pandas as pd
import os, statistics
import argparse
from scipy.signal import savgol_filter
from modbampredictnucsw import getThreshFromFile, predictOpenChromatin, predictNucleosomePos

parser = argparse.ArgumentParser(description='Visualizing modified bam files and nucleosome bed files',
                                 usage='python[3+] modbamvisualization.py -b file.bam OR -m file-moddata.txt -r "chrN:startpos-stoppos" [options]')

parser.add_argument('-b', '--bamfile', action='store', dest='b',
                    help='modified bam file, should be indexed. Can be filtered to locus or not. Optimal if name does not include dashes.')
parser.add_argument('-y', '--modtype', action='store', dest='y',
                    help='type of modification to visualize. This is only required if you are processing a .bam file. Allowable codes are: 5mC, 5hmC, 5fC, 5caC, 5hmU, 5fU, 5caU, 6mA, 8oxoG, Xao')
parser.add_argument('-m', '--moddatafile', action='store', dest='m',
                    help='data file from earlier run of this program, contains readnames and modification code for each position')
parser.add_argument('-n', '--nucbed', action='store', dest='n',
                    help='nucleosomes.bed file from cawlr sma or from modbampredictnuc-sw.py')
parser.add_argument('-k', '--pca', action='store_true', dest='k',
                    help='whether to output readnames to clusters file and plot PCA and kmeans clustering (the plot is only if using kmeans clustering in -c)')
parser.add_argument('-p', '--plot', action='store_true', dest='p',
                    help='whether to plot reads with modified positions, will cluster, outputs pdf')
parser.add_argument('-o', '--overlay', action='store_true', dest='o',
                    help='whether to plot reads with modifications overlayed with nucleosome positions, required -n')
parser.add_argument('-x', '--remove', action='store_true', dest='x',
                    help='whether to remove highly/chaotically modified reads. This will only work if you specify a threshold.')
parser.add_argument('-t', '--threshold', action='store', dest='t',
                    help='threshold between 0 and 1 to binarize modifications, mod score above this are called as true. Reccommend running predictthreshold.py or looking at score dist of positive and negative controls to determine this. If this is not set, will not binarize for plotting and the plotting will run much slower.')
parser.add_argument('-c', '--clusters', action='store', dest='c', default='2',
                    help='Options: manual, promoter, or any number n, "manual" means you will provide a file with premade cluster assignments using -e, "promoter" means to do manual clustering by promoter locations, requires -g. A number means do kmeans clustering with this many clusters, default 2, reccommend looking at PCA to help determine')
parser.add_argument('-r', '--region', action='store', dest='r',
                    help='region to calculate modification for and plot "chrN:startpos-stoppos"')
parser.add_argument('-g', '--genes', action='store', dest='g',
                    help='gene annotation in .bed format for visualization; optional. Bed file format should be one line per gene; fields beyond the sixth column will not be processed')
parser.add_argument('-a', '--aggregatenucpred', action='store', dest='a',
                    help='plot one prediction of nucleosome positions for each cluster. Must provide threshold file generated from predictthreshold.py containing predicted threshold, low threshold, and pos prob. You can make this file yourself as long as you follow the same formatting.')
parser.add_argument('-s', '--smoothingtype', action='store', default='some', dest='s',
                    help='Allowed options: none, some, full: smoothing means modifications are visualized as wider blocks where possible. This makes them easier to see. No smoothing means you will see the modified positions exactly as represented in the bam file. This may slow plotting down a little')
parser.add_argument('-d', '--bedgraph', action='store', dest='d',
                    help='bedgraph file with nucleosome positions')
parser.add_argument('-e', '--clusterfile', action='store', dest='e',
                    help='tsv file where each line has a readname column and a cluster label column. Requires -c manual')

args = parser.parse_args()
#k make PCA/cluster
#p plot just mod positions
#o plot mod positions with nucleosomes
#x remove highly/chaotically modified reads
#t threshold (0-255) if specified, this will binarize your modification data for plotting
#r region: "chrN:startpos-stoppos"
#c number of clusters

print('done loading modules and parsing args')


typesOfMods = {'5mC':[('C', 0, 'm')], '5hmC': [('C', 0, 'h')], '5fC': [('C', 0, 'f')], '5caC': [('C', 0, 'c')], 
               '5hmU': [('T', 0, 'g')], '5fU': [('T', 0, 'e')], '5caU': [('T', 0, 'b')], 
               '6mA': [('A', 0, 'a'), ('A', 0, 'Y')], '8oxoG': [('G', 0, 'o')], 'Xao': [('N', 0, 'n')]}


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
numToc = {-1.:'grey', 0.:'blue', 1.:'orange', 2.:'green', 3.:'purple', 4.:'red', 5.:'yellow', 6.:'greenyellow', 7.:'aqua', 8.:'deeppink', 9.:'goldenrod', 10.:'darkorchid'}


if not args.r:
    raise Exception('must specify region')
if not (args.b or args.m):
    raise Exception('must provide .bam file or moddata.txt file')
# if args.o and not args.n:
#     raise Exception('if you want to overlay nucleosomes, you must provide a nucleosomes.bed file')

# numclust = int(args.c)
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

            readstdev.append((np.std(knownpos), sum(knownpos)/len(knownpos)))
            readnames.append(readname)
            readdir.append(dir)
            modmatrix.append(thismat)
        print('done loading modified positions')
elif args.b:
    if not os.path.isfile(args.b):
        raise Exception('bam file does not exist')
    else:
        fileprefix = '.'.join(args.b.split('.')[:-1]) + '-' + args.y
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

                    posstag = typesOfMods[args.y]
                    if s.is_reverse: posstag = [(x[0], 1, x[2]) for x in posstag]
                    ml = None
                    for t in posstag:
                        if t in s.modified_bases:
                            ml = s.modified_bases[t]
                            break
                    if not ml:
                        print(readname, 'does not have modification information', s.modified_bases.keys())
                        continue

                    if s.has_tag('MM'):
                        skippedBase = -1 if s.get_tag('MM').split(',', 2)[0][-1] == '?' else 0
                    elif s.has_tag('Mm'):
                        skippedBase = -1 if s.get_tag('Mm').split(',', 2)[0][-1] == '?' else 0
                    else:
                        continue

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
                    mlpos = set([x[0] for x in ml])
                    for i in seqApos:
                        if i not in mlpos:
                            ml.append((i, skippedBase))
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


genes = []
openAndAdjPos = []
if args.g: #and (args.p or args.o):
    for line in open(args.g):
        line = line.rstrip().split('\t')
        if line[0].lstrip('chr') == chr.lstrip('chr'):
            geneleft, generight = int(line[1]), int(line[2])
            if geneleft < qstart and qstart < generight < qend: geneleft = qstart
            elif qstart < geneleft < qend  and generight > qend: generight = qend
            elif geneleft < qstart and generight > qend: geneleft, generight = qstart, qend
            # if line[3] == 'PHO5': print(geneleft, generight, qstart, qend)
            if geneleft >= qstart and generight <= qend:
                # print(geneleft, generight, line[3])
                mydircode = 2 if line[5] == '+' else 1
                genes.append([geneleft, generight-geneleft, mydircode, line[3]])
    print('done loading genes')
    # print(genes)
    genes.sort()

    if args.c == 'promoter': ###identify promoter regions to cluster on
        for i in range(len(genes)):
            if genes[i][2] == 2 and genes[i][0] > qstart: #gene in positive direction
                promend = genes[i][0] - qstart
                promstart = min((qstart, genes[i][0] - 1000))-qstart if i == 0 else max((genes[i][0] - min((1000, genes[i][0] - (genes[i-1][0]+genes[i-1][1]))), qstart)) - qstart
                genecompstart = genes[i][0] - qstart
                genecompend = min(genes[i][0] + min((1000, genes[i][1])), qend) - qstart
                openAndAdjPos.append(((promstart, promend), (genecompstart, genecompend)))
            elif genes[i][2] == 1 and genes[i][0] + genes[i][1] < qend: #gene is in negative direction
                promstart = genes[i][0] + genes[i][1] - qstart
                promend = min((qend, genes[i][0] + genes[i][1] + 1000)) -qstart if i == len(genes) - 1 else min((genes[i+1][0], genes[i][0] + genes[i][1] + 1000, qend)) - qstart
                genecompend = genes[i][0] + genes[i][1] - qstart
                genecompstart = max((qstart, genes[i][0] + genes[i][1] -1000)) - qstart if i ==0 else max((genes[i][0] + genes[i][1] -1000, genes[i-1][0]+genes[i-1][1], qstart)) - qstart    # genes[i][0] + min((1000, genes[i][1]))
                openAndAdjPos.append(((promstart, promend), (genecompstart, genecompend)))



clusterlabels = []
clusterReads = {}
clustertoreads = {}


if args.k or args.p or args.o or args.a: ##getting clustering done
    modmatrix = np.array(modmatrix)

    if args.x or args.c.isdigit(): ###if want to remove highly methylated reads
        pca = PCA(n_components=2)
        pca_features = pca.fit_transform(modmatrix)
        pca_df = pd.DataFrame(data=pca_features, columns=['PC1', 'PC2'])

        kmeans = KMeans(init="k-means++", n_clusters=2, n_init=4)
        kmeans.fit(pca_features)
        clusterlabels = kmeans.labels_.astype(float)

        subclust = [[] for i in range(2)]
        subclustreadnames = [[] for i in range(2)]
        readnametopos = {}
        for i in range(len(modmatrix)):
            subclust[int(clusterlabels[i])].append(modmatrix[i])
            subclustreadnames[int(clusterlabels[i])].append(readnames[i])
            readnametopos[readnames[i]] = i

        clusterlabels = [-1 for i in range(len(readnames))]
        fig, axs = plt.subplots(1, 2)
        for i in range(2):
            ###cluster based on methylation levels
            thismat = np.array(subclust[i])
            pca_features = pca.fit_transform(thismat)
            pca_df = pd.DataFrame(data=pca_features, columns=['PC1', 'PC2'])
            kmeans.fit(pca_features)
            theselabels = kmeans.labels_.astype(float)

            ###identify highly methylated cluster
            c1 = [sum(subclust[i][x]) for x in range(len(theselabels)) if theselabels[x] == 0.0]
            c2 = [sum(subclust[i][x]) for x in range(len(theselabels)) if theselabels[x] == 1.0]
            c1avg, c2avg = sum(c1)/len(c1), sum(c2)/len(c2)
            goodlabel = 0.0 if c1avg < c2avg else 1.0

            ####Remove highly methylated cluster
            goodnames = [subclustreadnames[i][x] for x in range(len(theselabels)) if theselabels[x] == goodlabel]
            if args.c == 'promoter' or args.c == 'manual':
                for j in range(len(goodnames)):
                    clusterlabels[readnametopos[goodnames[j]]] = 1
            elif args.c.isdigit():
                numclust = int(args.c)
                kmeans2 = KMeans(init="k-means++", n_clusters=numclust, n_init=4)
                goodmat = [subclust[i][x] for x in range(len(theselabels)) if theselabels[x] == goodlabel]
                ####Better cluster the remaining cluster based on user-defined number of clusters
                thismat = np.array(goodmat)
                pca_features = pca.fit_transform(thismat)
                pca_df = pd.DataFrame(data=pca_features, columns=['PC1', 'PC2'])
                kmeans2.fit(pca_features)
                theselabels = kmeans2.labels_.astype(float)
                for j in range(len(goodnames)):
                    thislabel = str(int(i)) + '.' + str(int(theselabels[j]))
                    clusterlabels[readnametopos[goodnames[j]]] = thislabel
                    if args.a:  # want to predict nucleosomes for clusters
                        if thislabel not in clusterReads: clusterReads[thislabel] = [0 for x in range(plotrange)]
                        clusterReads[thislabel] = [x + y for x, y in zip(clusterReads[thislabel], modmatrix[readnametopos[goodnames[j]]])]
                axs[i].scatter(pca_features[:, 0], pca_features[:, 1], c=[numToc[x] for x in theselabels], alpha=0.3)

        if args.k and not args.e:
            plt.savefig('-'.join([fileprefix, chr, str(qstart), str(qend)]) + '-kmeans-cluster.png', dpi=600)

        ####Remove reads in low clusters from matrix
        modmatrix = [modmatrix[x] for x in range(len(clusterlabels)) if clusterlabels[x] != -1]
        readdir = [readdir[x] for x in range(len(clusterlabels)) if clusterlabels[x] != -1]
        readnames = [readnames[x] for x in range(len(clusterlabels)) if clusterlabels[x] != -1]
        if args.c.isdigit(): clusterlabels = [x for x in clusterlabels if x != -1]
        elif args.c == 'promoter' or args.c == 'manual': clusterlabels = []

    if args.c == 'manual':
        readtoclusterlabel = {}
        for line in open(args.e):
            r, cl = line.rstrip().split('\t')
            readtoclusterlabel[r] = cl
        for i in range(len(readnames)):
            r = readnames[i]
            thislabel = readtoclusterlabel[r]
            clusterlabels.append(thislabel)
            if args.a:  # want to predict nucleosomes for clusters
                if thislabel not in clusterReads: clusterReads[thislabel] = [0 for x in range(plotrange)]
                clusterReads[thislabel] = [x + y for x, y in zip(clusterReads[thislabel], modmatrix[i])]

    elif args.c == 'promoter':
        modchunks, adjchunks = [], []
        for i in range(len(modmatrix)):
            thismat = modmatrix[i]
            regiongroups = []
            currstart = 0
            for region in openAndAdjPos:
                modinprom = []
                currstart = region[0][0]
                while currstart < region[0][1]-25:
                    modprob = thismat[currstart:currstart+25]
                    onlyAs = [x for x in modprob if x != -1]
                    # if currstart % 100 == 0: modpeaks.append(sum(onlyAs))
                    if len(onlyAs) >= 5:
                        if sum(onlyAs) / len(onlyAs) >= threshold: modinprom.append(1)
                        else: modinprom.append(0)
                    currstart += 1

                if len(modinprom) == 0:
                    regiongroups.append('0')
                    continue

                bestchunks = []
                for j in range(0, len(modinprom) - 500, 100):
                    modprob = modinprom[j:j + 500]
                    bestchunks.append(sum(modprob)/len(modprob))
                modprob = modinprom[max(0,len(modinprom)-500):len(modinprom)]
                bestchunks.append(sum(modprob)/len(modprob))

                chunksthreshold = 0.1 if threshold > 0.7*255 else 0.02
                promopencode = 1 if max(bestchunks) >= chunksthreshold else 0
                ###This is for R10 max(bestchunks) >= 0.02 else 0
                ###This is for R9 max(bestchunks) >= 0.1 else 0
                # print(promopencode, bestchunks)
                regiongroups.append(str(promopencode))
            thislabel = '.'.join(regiongroups)
            clusterlabels.append(thislabel)
            if args.a:  # want to predict nucleosomes for clusters
                if thislabel not in clusterReads: clusterReads[thislabel] = [0 for x in range(plotrange)]
                clusterReads[thislabel] = [x + y for x, y in zip(clusterReads[thislabel], thismat)]
            if '.'.join(regiongroups) not in clustertoreads: clustertoreads['.'.join(regiongroups)] = 0
            clustertoreads['.'.join(regiongroups)] += 1
        for i in clustertoreads:
            print(i, clustertoreads[i])

    if args.k and not args.e:
        clusterout = open(fileprefix + '-readnamesToClusters.tsv', 'w')
        for i in range(len(readnames)):
            clusterout.write(readnames[i] + '\t' + str(clusterlabels[i]) + '\n')
        clusterout.close()
    print('done clustering')





clusterToNucPos = {}
combinedNucPred = []
if args.a:
    clusterReads = [[key, item] for key, item in clusterReads.items()]
    threshold, lowthresh, posprob = getThreshFromFile(args.a)
    clusterReads.sort(key=lambda x:sum(x[1]))
    bestclust = [0 for x in range(plotrange)]
    for i in range(len(clusterReads)):
        bestclust = [x+y for x, y in zip(bestclust, clusterReads[i][1])]

        clusterReads[i][1] = [-1 if x < 0 else math.log2((x/40)+1) for x in clusterReads[i][1]]

        clustermax = max(clusterReads[i][1])

        clusterReads[i][1] = [-1 if x < 0 else int((x / clustermax) * 255) for x in clusterReads[i][1]]

    bestclust = [-1 if x < 0 else math.log2((x / 1000) + 1) for x in bestclust]
    bestclustmax = max(bestclust)
    bestclust = [-1 if x < 0 else int((x / bestclustmax) * 255) for x in bestclust]

    # bestclust = [-1 if x < 0 else int((x**2 / bestclustmax) * 255) for x in bestclust]
    # bestclust = [-1 if x < 0 else (0.0626 * x) ** 2 for x in bestclust]
    # bestclustmax = max(bestclust)
    # bestclust = [-1 if x < 0 else int((x / bestclustmax) * 255) for x in bestclust]

    thisposprob = posprob - 0.1 if posprob - .1 > 0.5 else posprob + 0.1
    bestprednuc = predictOpenChromatin([bestclust], plotrange, 25, lowthresh, threshold, thisposprob, False)
    finalbestpred = predictNucleosomePos(bestprednuc, plotrange, True)
    thisposprob = 0.5 if posprob < 0.5 else posprob
    predictednuc = predictOpenChromatin([x[1] for x in clusterReads], plotrange, 25, lowthresh, threshold, thisposprob, False)
    finalnucpred = predictNucleosomePos(predictednuc, plotrange, True)

    for i in range(len(clusterReads)):
        clusterToNucPos[clusterReads[i][0]] = finalnucpred[i]
    combinedNucPred = finalbestpred[0]
    print('done predicting nucleosome positions for clusters')

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
        nomodcolor = 1 if readdir[z] == '-' else 2
        if args.s == 'some': #smoothing allowed
            currColor = 0
            lasta, lastpos = 0, 0
            colbybp = [0 for x in range(plotrange)]
            for i in range(qstart, qend):
                modprob = thesemods[i - qstart]
                if modprob != -1:
                    if i - lasta > 5:
                        thesebounds = [i - 3, i + 3]
                    elif i - lasta > 3:
                        thesebounds = [int(i - ((i - lasta) / 2)), i + 3]
                    else:
                        thesebounds = [(i - (i - lasta)) + 1, i + 3]

                    lasta = i
                    if args.t:
                        if modprob >= threshold:
                            currColor = 3  # red
                        else:
                            currColor = nomodcolor
                    else:
                        currColor = modprob
                else: ###never smooth whitespace
                    thesebounds = [i, i+1]
                    currColor = 0
                for j in range(thesebounds[0], thesebounds[1]):
                    if j < qend:
                        colbybp[j - qstart] = currColor  # codetocol[currColor]
        elif args.t: #thresholded
            if args.s == 'none':
                colbybp = [3 if thesemods[x] >= threshold else nomodcolor for x in range(plotrange)]
            ###BELOW IS ACTUALLY ALTERNATE SMOOTHING METHOD
            elif args.s == 'full':
                colbybp = [0 for x in range(plotrange)]
                i = 25
                currstart = 0
                while i < plotrange:
                    modprob = thesemods[currstart:i]
                    onlyAs = [x for x in modprob if x != -1]
                    if len(onlyAs) >= 5:
                        if sum(onlyAs)/len(onlyAs) >= threshold: currColor = 3
                        else: currColor = nomodcolor
                        for j in range(currstart, i):
                            if colbybp[j] != 3: colbybp[j] = currColor
                    i += 1
                    currstart += 1


        if args.t: ###removing small bits of whitespace
            lastcol, laststart = 0, 0
            for i in range(plotrange):
                if colbybp[i] != lastcol:
                    if lastcol == 0 and i - laststart <= 5:
                        for j in range(laststart, i):
                            colbybp[j] = nomodcolor
                    lastcol, laststart = colbybp[i], i

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
        if thismax > 0:
            clusteravg[label] = [x/thismax for x in clusteravg[label]]
    readcolbybp.sort()#(key=lambda x: [sum(clusteravg[x[0]]), -1* sum(x[1])])
    print('done calculating cluster averages and smoothing if desired')


nucfrombedgraph = [0 for x in range(plotrange)]
if args.d:
    for line in open(args.d):
        if line[:5] != 'track':
            line = line.rstrip().split()
            if line[0] == chr:
                if qstart <= int(line[2]) < qend:
                    nucfrombedgraph[int(line[2])-qstart] = float(line[3])
    thismax = max(nucfrombedgraph)
    nucfrombedgraph = [x/thismax for x in nucfrombedgraph]
    # print(nucfrombedgraph)
    print('done processing bedgraph file')


if args.p: #plot mod pos
    totreads = len(readcolbybp) + (len(clusteravg) * 2) + len(genes)  # num clusters
    if args.d: totreads += 2
    figheight = 100 if totreads/5 > 98 else 2+int(totreads/5)
    plt.figure(figsize = (12, figheight))
    ax = plt.axes(frameon=False)
    box = ax.get_position()
    ax.set_position([box.x0 + (box.width*0.22), box.y0, box.width * 0.85, box.height])
    ticklocs, ticklabels = [], []
    c = 0
    xaxis = [x+qstart for x in range(plotrange)]

    for g in genes:
        rectangle1 = mplpatches.Rectangle((g[0], totreads-c), g[1], 0.9, facecolor=codetocol[g[2]],
                                          linewidth=0, alpha=0.9)
        ax.add_patch(rectangle1)
        ticklocs.append(0.5 + totreads - c)
        ticklabels.append(g[3])
        if g[2] == 2 and g[0] > qstart:       #gene in positive direction
            plt.arrow(max((qstart, g[0]-int(plotrange/20))), 0.5+totreads-c, int(plotrange/20), 0, length_includes_head=True, head_width=0.3, head_length=int(plotrange/100), facecolor='black')
        elif g[2] == 1 and g[0] + g[1] < qend:
            plt.arrow(min((qend, g[0]+g[1]+int(plotrange/20))), 0.5+totreads-c, -1*int(plotrange/20), 0, length_includes_head=True, head_width=0.3, head_length=int(plotrange/100), facecolor='black')
        c += 1

    if args.d:
        c += 1
        plt.plot(xaxis, [(totreads-c)+(x * 1.6) + 0.2 for x in nucfrombedgraph], color='black')
        ticklocs.append(0.5 + totreads - c)
        ticklabels.append(args.d.split('/')[-1].split('.')[0])
        c += 1

    lastcluster = -2
    for colbybp in readcolbybp:
        if colbybp[0] != lastcluster:
            c += 1
            plt.plot(xaxis, [(totreads-c)+(x*1.6) +0.2 for x in clusteravg[colbybp[0]]], color='black')
            ticklocs.append(0.5 + totreads - c)
            ticklabels.append('Cluster ' + str(colbybp[0]))
            c += 1
            lastcluster = colbybp[0]
        mydircode = 2 if colbybp[3] == '+' else 1
        readname = colbybp[2] + ' ' + colbybp[3] #+ ' ' + ' '.join(colbybp[4]) #str(round(colbybp[4], 3))
        colbybp = colbybp[1]
        rectshapes = []
        laststart, lastcol = 0, 0
        for i in range(plotrange):
            if colbybp[i] != lastcol or i == plotrange-1:
                rectshapes.append([laststart+qstart, i-laststart, lastcol, totreads-c])
                laststart, lastcol = i, colbybp[i]
        ticklocs.append(0.5 + totreads - c)
        ticklabels.append(readname)
        c += 1

        lastcol = 0
        for i in range(len(rectshapes)):
            rect = rectshapes[i]
            if args.t:
                if rect[1] <= 5 and rect[2] == 0: rect[2] = mydircode
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
    plotname.append(args.s + 'Smoothing')
    plt.savefig('-'.join(plotname) + '-modifiedReads.pdf', format='pdf', dpi=600)
    print('done plotting modified reads')





if args.o: #plot mod pos with overlaid nucleosomes
    totreads = len(readcolbybp) + (len(clusteravg) * 2) + len(genes)  # num clusters
    if args.a: totreads += len(clusteravg) * 2 + 1
    if args.d: totreads += 2
    figheight = 100 if int(totreads/5) > 98 else int(totreads/5)+2
    plt.figure(figsize = (12, figheight))
    ax = plt.axes(frameon=False)
    box = ax.get_position()
    ax.set_position([box.x0 + (box.width * 0.22), box.y0, box.width * 0.85, box.height])
    ticklocs, ticklabels = [], []
    c = 0
    xaxis = [x+qstart for x in range(plotrange)]

    for g in genes:
        rectangle1 = mplpatches.Rectangle((g[0], totreads - c), g[1], 0.9, facecolor=codetocol[g[2]],
                                          linewidth=0, alpha=0.9)
        ax.add_patch(rectangle1)
        ticklocs.append(0.5 + totreads - c)
        ticklabels.append(g[3])
        if g[2] == 2 and g[0] > qstart:  # gene in positive direction
            plt.arrow(max((qstart, g[0] - int(plotrange / 20))), 0.5 + totreads - c, int(plotrange / 20), 0,
                      length_includes_head=True, head_width=0.3, head_length=int(plotrange / 100), facecolor='black')
        elif g[2] == 1 and g[0] + g[1] < qend:
            plt.arrow(min((qend, g[0] + g[1] + int(plotrange / 20))), 0.5 + totreads - c, -1 * int(plotrange / 20), 0,
                      length_includes_head=True, head_width=0.3, head_length=int(plotrange / 100), facecolor='black')
        c += 1

    if args.d:
        c += 1
        plt.plot(xaxis, [(totreads - c) + (x * 1.6) + 0.2 for x in nucfrombedgraph], color='black')
        ticklocs.append(0.5 + totreads - c)
        ticklabels.append(args.d.split('/')[-1].split('.')[0])
        c += 1

    if args.a:
        for block in combinedNucPred:
            rectangle1 = mplpatches.Rectangle((block[0] + qstart, totreads - c), block[1], 0.9,
                                              facecolor='grey', linewidth=0, alpha=0.75)
            ax.add_patch(rectangle1)
        ticklocs.append(0.5 + totreads - c)
        ticklabels.append('Combined nucleosome predictions')
        c += 1
    lastcluster = -2
    for colbybp in readcolbybp:
        if colbybp[0] != lastcluster:
            c += 1
            plt.plot(xaxis, [(totreads-c)+(x*1.6)+0.2 for x in clusteravg[colbybp[0]]], color='black')
            ticklocs.append(0.5+ totreads-c)
            ticklabels.append('Cluster ' + str(colbybp[0]))
            c += 1
            if args.a:
                for block in clusterToNucPos[colbybp[0]]:
                    rectangle1 = mplpatches.Rectangle((block[0]+qstart, totreads - c + .3), block[1], 0.6,
                                                      facecolor='grey', linewidth=0, alpha=0.75)
                    ax.add_patch(rectangle1)
                    ticklocs.append(0.5 + totreads - c)
                    ticklabels.append('Nucleosome predictions for cluster ' + str(colbybp[0]))
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
        if args.n:
            for block in readToNucPos[readname.split(' ')[0]]:
                rectangle1 = mplpatches.Rectangle((block[0], totreads-c), block[1]-block[0], 0.9, facecolor='grey',
                                                  linewidth=0, alpha=0.5)
                ax.add_patch(rectangle1)
        ticklocs.append(0.5 + totreads-c)
        ticklabels.append(readname)
        c += 1

        # if c == 3:
        #     print(colbybp)
        #     print(rectshapes)

        for i in range(len(rectshapes)):
            rect = rectshapes[i]
            if args.t:
                if rect[1] <= 5 and rect[2] == 0: rect[2] = mydircode
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
    plotname.append('modifiedReads')
    if args.n: plotname.append(args.n.split('-')[-1].split('.')[0])
    else: plotname.append('clusternucpred')
    plotname.append(args.s+'Smoothing')
    plt.savefig('-'.join(plotname)  + '.pdf', format='pdf', dpi=600)
    print('done plotting modified reads overlaid with nucleosomes')
