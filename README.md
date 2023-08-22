# modbamvisualization
Visualizing modifications from .bam file



Visualizing modified bam files and nucleosome bed files, as well as doing some simple nucleosome prediction

## Table of Contents

- Options
- Visualization basics
- Understanding modification score distributions and thresholds
- Predicting and visualizing nucleosomes
  - Nucleosome prediction options
  - Visualizing nucleosomes
- Aggregating modified position calling as a bedgraph
  - Calling peaks from this data

## Options
```bash
usage: python[3+] modbamvisualization.py -b file.bam OR -m file-moddata.txt -r "chrN:startpos-stoppos" [options]

Visualizing modified bam files and nucleosome bed files

options:
  -h, --help            show this help message and exit
  -b B, --bamfile B     modified bam file, should be indexed. Can be filtered to locus or not.
                        Optimal if name does not include dashes.
  -y Y, --modtype Y     type of modification to visualize. This is only required
                        if you are processing a .bam file. Allowable codes
                        are: 5mC, 5hmC, 5fC, 5caC, 5hmU, 5fU, 5caU, 6mA, 8oxoG, Xao
  -m M, --moddatafile M
                        data file from earlier run of this program,
                        contains readnames and modification code for each position
  -n N, --nucbed N      nucleosomes.bed file from cawlr sma or from modbampredictnuc-sw.py
  -k, --pca             whether to output readnames to clusters file and plot PCA and
                        kmeans clustering (the plot is only if using kmeans clustering in -c)
  -p, --plot            whether to plot reads with modified positions, will cluster, outputs pdf
  -o, --overlay         whether to plot reads with modifications overlayed
                        with nucleosome positions, required -n
  -x, --remove          whether to remove highly/chaotically modified reads.
                        This will only work if you specify a threshold.
  -t T, --threshold T   threshold between 0 and 1 to binarize modifications,
                        mod score above this are called as true. Reccommend
                        running predictthreshold.py or looking at score dist
                        of positive and negative controls to determine this. If
                        this is not set, will not binarize for plotting and the plotting will run much slower.
  -c C, --clusters C    Options: manual, promoter, or any number n,
                        "manual" means you will provide a file with premade cluster assignments using -e,
                        "promoter" means to do manual clustering by promoter locations, requires -g.
                        A number means do kmeans clustering with this many clusters, default 2, reccommend looking at PCA to help determine
  -r R, --region R      region to calculate modification for and plot "chrN:startpos-stoppos"
  -g G, --genes G       gene annotation in .bed format for visualization; optional.
                        Bed file format should be one line per gene;
                        fields beyond the sixth column will not be processed
  -a A, --aggregatenucpred A
                        plot one prediction of nucleosome positions for each cluster.
                        Must provide threshold file generated from
                        predictthreshold.py containing predicted threshold, low threshold, and pos prob.
                        You can make this file yourself as long as you follow the same formatting.
  -s S, --smoothingtype S
                        Allowed options: none, some, full:
                        smoothing means modifications are visualized as wider blocks where possible.
                        This makes them easier to see. No smoothing means you will see the modified
                        positions exactly as represented in the bam file. This may slow plotting down a little
  -d D, --bedgraph D    bedgraph file with nucleosome positions
  -e E, --clusterfile E
                        tsv file where each line has a readname column and a cluster label column. Requires -c manual
```

## Visualization basics

I reccommend starting by running something like:

```bash
python modbamvisualization.py -b file.bam -r "chrN:startpos-stoppos" -k
```

this will generate two files: file-moddata.txt and file-kmeans-cluster.png

You can look at the clustered PCA to decide how many clusters you want to move forward with.

For example, in the below image, each panel is a strand, and you can see that each strand clusters in about 3 clusters, likely due to different densities of modification.

![alt text](https://github.com/cafelton/modbamvisualization/blob/main/r9_nuc_mod_region_pho5-chrII-430000-435000-kmeans-cluster-subcluster.png?raw=true)

When you've looked at that and made sure everything looks ok, you can rerun again with:

```bash
python modbamvisualization.py -m file-moddata.txt -r "chrN:startpos-stoppos" -p -c 4
```

This will generate a plot like below, with clustered reads and modifications in red. The strands are slightly different shades of blue. Running like this is much slower and does not produce the clearest image though. I actually reccommend skipping straight to the next steps.

![alt text](https://github.com/cafelton/modbamvisualization/blob/main/no-threshold-reads.png?raw=true)

## Understanding modification score distributions and thresholds

If you want a clearer image, I reccommend setting a threshold. 0.7 tends to be a good starting place. If you have positive and negative control data, you can run predictthreshold.py as below. This will output a predictedThreshold.tsv data file, the first line of which is a suggested threshold to move forward with for visualization.

```bash
python predictthreshold.py -p pos.bam -n neg.bam
```

Note that the third line of the output file is a good shorthand for how good the separation between your positive and negative controls is - the closer this number is to 1, the better distinction there is between open and closed chromatin and the better the nucleosome predictions will be.

You can also run this without control data and instead use -s sample.bam, but this is much less robust and not reccommended.

If you want to visualize the score distribution of any number of files, you can run:

```bash
python makeScoreDistFromModBam.py pos.bam neg.bam [etc.bam]
```

![alt text](https://github.com/cafelton/modbamvisualization/blob/main/ScoreDist.png?raw=true)

Once you figure out a threshold, you can visualize the data with it and run:

```bash
python modbamvisualization.py -m file-moddata.txt -r "chrN:startpos-stoppos" -p -c 4 -t 0.7
```

You can see that the clusters have clearer peaks in open chromatin now.

![alt text](https://github.com/cafelton/modbamvisualization/blob/main/threshold-reads.png?raw=true)

## Predicting and visualizing nucleosomes

If you already have a nucleosome prediction bed file generatd by cawlr sma or another tool, or any other bed file with the read name corresponding with the read names in your bam file that you wish to visualize, you can skip to the end of this section.

If you need to predict nucleosomes, first run predictthreshold.py as described above or manually make a .tsv file in the correct format. Below is an example if you don't know where to start:

```bash
predicted threshold:	0.8
lowest possible threshold:	0.7
pos nuc prob above threshold:	0.66
```

You can then run the nucleosome prediction with a command like below:

```bash
python modbampredictnuc-sw.py -m file-moddata.txt -t predictedThreshold.tsv -r chr:startpos-stoppos
```

### Nucleosome prediction options

```bash
usage: python[3+] modbamvisualization.py -m file-moddata.txt -r "chrN:startpos-stoppos" -t predictedThreshold.tsv [options]

Predicting nucleosome positions using a sliding window

options:
  -h, --help            show this help message and exit
  -m M, --moddatafile M
                        data file from modbedvisualization, contains readnames and modification code for each position. The region this
                        was run with should be the same as the inputted region.
  -t T, --thresholdfile T
                        file generated from predictthreshold.py containing predicted threshold, low threshold, and pos prob. You can
                        make this file yourself as long as you follow the same formatting.
  -r R, --region R      region to calculate nucleosomes for. make sure this matches the region modbedvisualization was run with.
                        "chrN:startpos-stoppos"
  -w W, --windowsize W  window size to scan, larger is more stringent for searching for larger regions of open chromatin. Default: 30,
                        unit is base pairs
  -o, --predopenonly    whether to only predict open chromatin and not nucleosomes, default: False
```

### Visualizing nucleosomes

To actually visualize nucleosomes, you just add the -o and -n options to your visualization command as below:

```bash
python modbamvisualization.py -m file-moddata.txt -r "chrN:startpos-stoppos" -c 3 -t 0.7 -o -n file_nucleosomes.bed
```

The open chromatin can be seen where there are red modifications and a lack of grey nucleosomes.

![alt text](https://github.com/cafelton/modbamvisualization/blob/main/nuc-threshold.png?raw=true)

If you also want to generate nucleosome predictions for each cluster and for everything combined, you can do that within modbamvisualization by running:

```bash
python modbamvisualization.py -m file-moddata.txt -r "chrN:startpos-stoppos" -c 2 -t 0.7 -o -n file_nucleosomes.bed -a predictedThreshold.tsv
```

You can also run the -a option without having already predicted nucleosome positions for the individual reads:

```bash
python modbamvisualization.py -m file-moddata.txt -r "chrN:startpos-stoppos" -c 2 -t 0.7 -o -a predictedThreshold.tsv -g genes.bed
```
![alt text](https://github.com/cafelton/modbamvisualization/blob/main/plotwithclusternucpred.png?raw=true)

## Aggregating modified position calling as a bedgraph

```bash
usage: python[3+] modbamvisualization.py -b file.bam OR -m file-moddata.txt -r "chrN:startpos-stoppos" [options]

Visualizing modified bam files and nucleosome bed files

options:
  -h, --help            show this help message and exit
  -b B, --bamfile B     modified bam file, should be indexed. Can be filtered to locus or not. Optimal if name does not include
                        dashes.
  -y Y, --modtype Y     type of modification to visualize. This is only required if you are processing a .bam file. Allowable
                        codes are: 5mC, 5hmC, 5fC, 5caC, 5hmU, 5fU, 5caU, 6mA, 8oxoG, Xao
  -r R, --region R      region to calculate modification for and plot "chrN:chrlen" OR "chrN:startpos-stoppos"
  -t T, --threshold T   threshold between 0 and 1 to binarize modifications, mod score above this are called as true.
                        Reccommend running predictthreshold.py or looking at score dist of positive and negative controls to
                        determine this.
  -e E, --clusterfile E
                        tsv file where each line has a readname column and a cluster label column. Requires -c manual
  -s, --smoothing       whether to smooth modifications curve
```

### Calling peaks from this data

```bash
usage: python[3+] modbamvisualization.py -b file.bam OR -m file-moddata.txt -r "chrN:startpos-stoppos" [options]

Visualizing modified bam files and nucleosome bed files

options:
  -h, --help           show this help message and exit
  -b B, --bedgraph B   bedgraph file
  -t T, --threshold T  threshold between 0 and 1 to binarize modifications, mod score above this are called as true. Reccommend
                       running predictthreshold.py or looking at score dist of positive and negative controls to determine
                       this.
  -r R, --region R     region to calculate modification for and plot "chrN:chrlen" OR "chrN:startpos-stoppos"
```
