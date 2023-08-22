import sys
import pysam
import argparse



parser = argparse.ArgumentParser(description='Visualizing modified bam files and nucleosome bed files',
                                 usage='python[3+] modbamvisualization.py -b file.bam OR -m file-moddata.txt -r "chrN:startpos-stoppos" [options]')

parser.add_argument('-b', '--bedgraph', action='store', dest='b',
                    help='bedgraph file')
parser.add_argument('-t', '--threshold', action='store', dest='t',
                    help='threshold between 0 and 1 to binarize modifications, mod score above this are called as true. Reccommend running predictthreshold.py or looking at score dist of positive and negative controls to determine this.')
parser.add_argument('-r', '--region', action='store', dest='r',
                    help='region to calculate modification for and plot "chrN:chrlen" OR "chrN:startpos-stoppos"')

args = parser.parse_args()


if '-' not in args.r:
    chr, plotrange = args.r.split(':')
    plotrange = int(plotrange)
    qstart = 0
else:
    chr, temp = args.r.split(':')
    qstart, qend = temp.split('-')
    qstart, qend = int(qstart), int(qend)
    plotrange = qend - qstart

threshold = int(float(args.t)*255)
print(threshold)

fileprefix = '.'.join(args.b.split('.')[:-1])

modpos = [-1 for x in range(plotrange)]

for line in open(args.b):
    line = line.rstrip().split('\t')
    if line[0] == chr:
        modpos[int(line[2])-qstart] = float(line[3])
print('done processing file')
calledregions = [0 for x in range(plotrange)]

for i in range(0, plotrange-50, 5):
    thismodpos = modpos[i:i + 50]
    onlyAs = [x for x in thismodpos if x != -1]
    # if 377906<i<377995:
    #     print('closed', sum(onlyAs) / len(onlyAs), onlyAs)
    # if 378900 < i < 379000:
    #     print('open', sum(onlyAs)/len(onlyAs), onlyAs)
    if len(onlyAs) >= 5:
        if sum(onlyAs) / len(onlyAs) >= threshold:
            for j in range(i, i+50):
                calledregions[j] = 1
print('done finding peaks')
lastcol, laststart = 0, 0
for i in range(plotrange):
    if calledregions[i] != lastcol:
        if lastcol == 0 and i - laststart <= 10:
            for j in range(laststart, i):
                calledregions[j] = 1
        lastcol, laststart = calledregions[i], i
print('done smoothing')
out = open(fileprefix + '-peaks.bed', 'w')
laststart, lastcol = 0, 0
for i in range(plotrange):
    if calledregions[i] != lastcol or i == plotrange-1:
        if lastcol == 1:
            thisregionAs = [x for x in modpos[laststart:i] if x != -1]
            out.write('\t'.join([chr, str(laststart+qstart), str(i+qstart), 'peakwidth=' + str(i-laststart) + ' peakheight=' + str(round(1-((sum(thisregionAs)/len(thisregionAs))/255), 3))]) + '\n')
        laststart, lastcol = i, calledregions[i]
out.close()