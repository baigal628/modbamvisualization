import sys
import pysam
import argparse
from scipy.signal import savgol_filter
import numpy as np


parser = argparse.ArgumentParser(description='Visualizing modified bam files and nucleosome bed files',
                                 usage='python[3+] modbamvisualization.py -b file.bam OR -m file-moddata.txt -r "chrN:startpos-stoppos" [options]')

parser.add_argument('-b', '--bamfile', action='store', dest='b',
                    help='modified bam file, should be indexed. Can be filtered to locus or not. Optimal if name does not include dashes.')
parser.add_argument('-y', '--modtype', action='store', dest='y',
                    help='type of modification to visualize. This is only required if you are processing a .bam file. Allowable codes are: 5mC, 5hmC, 5fC, 5caC, 5hmU, 5fU, 5caU, 6mA, 8oxoG, Xao')
parser.add_argument('-r', '--region', action='store', dest='r',
                    help='region to calculate modification for and plot "chrN:chrlen" OR "chrN:startpos-stoppos"')
parser.add_argument('-t', '--threshold', action='store', dest='t',
                    help='threshold between 0 and 1 to binarize modifications, mod score above this are called as true. Reccommend running predictthreshold.py or looking at score dist of positive and negative controls to determine this.')
parser.add_argument('-e', '--clusterfile', action='store', dest='e',
                    help='tsv file where each line has a readname column and a cluster label column. Requires -c manual')
parser.add_argument('-s', '--smoothing', action='store_true', dest='s',
                    help='whether to smooth modifications curve')
args = parser.parse_args()


typesOfMods = {'5mC':[('C', 0, 'm')], '5hmC': [('C', 0, 'h')], '5fC': [('C', 0, 'f')], '5caC': [('C', 0, 'c')],
               '5hmU': [('T', 0, 'g')], '5fU': [('T', 0, 'e')], '5caU': [('T', 0, 'b')],
               '6mA': [('A', 0, 'a'), ('A', 0, 'Y')], '8oxoG': [('G', 0, 'o')], 'Xao': [('N', 0, 'n')]}
compbase = {'A':'T', 'T':'A', 'C':'G', 'G':'C', 'N':'N'}
def getcomp(seq):
    newseq = []
    for base in seq: newseq.append(compbase[base])
    return ''.join(newseq)#newseq[::-1]


if '-' not in args.r:
    chr, plotrange = args.r.split(':')
    plotrange = int(plotrange)
    qstart = 0
else:
    chr, temp = args.r.split(':')
    qstart, qend = temp.split('-')
    qstart, qend = int(qstart), int(qend)
    plotrange = qend - qstart

if args.t: threshold = int(float(args.t)*255)

readtoclusterlabel = {}
moddict = {}
if args.e:
    for line in open(args.e):
        r, cl = line.rstrip().split('\t')
        readtoclusterlabel[r] = cl
        if cl not in moddict: moddict[cl] =  {x:[] for x in range(plotrange)}
if not args.e:
    moddict['0'] = {x:[] for x in range(plotrange)}


fileprefix = '.'.join(args.b.split('/')[-1].split('.')[:-1]) + '-' + args.y
# out = open('-'.join([fileprefix, chr, str(qstart), str(qend)]) + '-moddata.txt', 'w')





samfile = pysam.AlignmentFile(args.b, "rb")

myfetch = samfile.fetch(chr) if '-' not in args.r else samfile.fetch(chr, qstart, qend)

for s in myfetch:
    if not s.is_secondary:
        alignstart, alignend = s.reference_start, s.reference_end
        readname = s.query_name
        if not args.e or readname in readtoclusterlabel:
            print(readname)
            cigar = s.cigartuples

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

            skippedBase = None
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
            posdict = dict(posOnGenome)
            for i in posdict:
                if posdict[i] != -1:
                    j = i-alignstart if '-' in args.r else i
                    if 0 <= j < plotrange:
                        cl = readtoclusterlabel[readname] if args.e else '0'
                        if args.t:
                            ismod = 1 if posdict[i] > threshold else 0
                            moddict[cl][j].append(ismod)
                        else:
                            moddict[cl][j].append(posdict[i])
samfile.close()
print('done parsing region in bam file')



for clust in moddict:
    thisout = fileprefix + '-cluster' + clust
    if args.t: thisout += '-' + str(args.t) + 'threshold'
    if args.s: thisout += '-smoothed'
    out = open(thisout + '.bedgraph', 'w')
    if not args.s:
        for i in range(plotrange):
            # moddict[clust][i] = [x for x in moddict[clust][i] if x > 0]
            if len(moddict[clust][i]) > 0:
                if args.t:
                    out.write('\t'.join([chr, str(qstart + i - 1), str(qstart + i), str(int((255*(sum(moddict[clust][i]) / len(moddict[clust][i])))))]) + '\n')
                else:
                    out.write('\t'.join([chr, str(qstart + i-1), str(qstart + i), str(int(sum(moddict[clust][i])/len(moddict[clust][i])))]) + '\n')
    else:
        thisavg = [-1 for x in range(plotrange)]
        lastpos, lastval = 0, 0
        for i in range(plotrange):
            if len(moddict[clust][i]) > 0:
                thisval = sum(moddict[clust][i])/len(moddict[clust][i])
                if i-lastpos > 1:
                    smoothedval = np.linspace(lastval, thisval, num=1+i-lastpos)
                    for j in range(1, i-lastpos):
                        thisavg[lastpos+j] = smoothedval[j]
                thisavg[i] = thisval
                lastpos, lastval = i, thisval
        # print(clust, thisavg[:30])
        # thisavg = savgol_filter(thisavg, 150, 2)
        # print(clust, thisavg[:30])
        for i in range(plotrange):
            if args.t:
                out.write('\t'.join([chr, str(qstart + i - 1), str(qstart + i), str(int((255*thisavg[i])))]) + '\n')
            else:
                out.write('\t'.join([chr, str(qstart + i-1), str(qstart + i), str(int(thisavg[i]))]) + '\n')
    out.close()