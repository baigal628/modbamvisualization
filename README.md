# modbamvisualization
Visualizing modifications from .bam file

usage: python[3+] modbamvisualization.py -b file.bam OR -m file-moddata.txt -r "chrN:startpos-stoppos" [options]

Visualizing modified bam files and nucleosome bed files

options:

  -h, --help            show this help message and exit
  
  -b B, --bamfile B     modified bam file, should be indexed. Can be filtered to locus or not. Optimal if name does not include dashes.
  
  -m M, --moddatafile M
                        data file from earlier run of this program, contains readnames and modification code for each position
  
  -n N, --nucbed N      nucleosomes.bed file from cawlr sma
 
  -k, --pca             whether to plot PCA and run kmeans clustering (default clusters: 3, specify otherwise with -c)
  
  -p, --plot            whether to plot reads with modified positions, will cluster, outputs pdf
  
  -o, --overlay         whether to plot reads with modifications overlayed with nucleosome positions, required -n
  
  -x, --remove          whether to remove highly/chaotically modified reads. This will only work if you specify a threshold.
 
  -t T, --threshold T   threshold between 0 and 1 to binarize modifications, mod score above this are called as true. Reccommend
                        looking at score dist of positive and negative controls to determine this. If this is not set, will not
                        binarize for plotting.
  
  -c C, --clusters C    number of clusters for plotting, deault 3, reccommend looking at PCA to help determine
  
  -r R, --region R      region to calculate modification for and plot "chrN:startpos-stoppos"



I reccommend starting by running something like:

python[3+] modbamvisualization.py -b file.bam -r "chrN:startpos-stoppos" -k

this will generate two files: file-moddata.txt and file-kmeans-cluster.png

You can look at the clustered PCA to decide how many clusters you want to move forward with

For example, the below image shows clearly that we have + strand reads (left), - strand reads (right) and less modifications (bottom) and more modifications (top)

![alt text](https://github.com/cafelton/modbamvisualization/blob/main/clusters.png?raw=true)

When you've looked at that and made sure everything looks ok, you can rerun again with:

python[3+] modbamvisualization.py -m file-moddata.txt -r "chrN:startpos-stoppos" -p -c 4

This will generate a plot like below, with clustered reads and modifications in red. The strands are slightly different shades of blue.

![alt text](https://github.com/cafelton/modbamvisualization/blob/main/no-threshold-reads.png?raw=true)

If you want a clearer image, I reccommend setting a threshold:

python[3+] modbamvisualization.py -m file-moddata.txt -r "chrN:startpos-stoppos" -p -c 4 -t 0.7

You can see that the clusters have clearer peaks in open chromatin now.

![alt text](https://github.com/cafelton/modbamvisualization/blob/main/threshold-reads.png?raw=true)

If you have run cawlr or have a nucleosome position bed file (one line of nucleosome positions per read), you can overlay the nucleosomes with the reads:

python[3+] modbamvisualization.py -m file-moddata.txt -r "chrN:startpos-stoppos" -c 4 -t 0.7 -o -n file_nucleosomes.bed

The open chromatin can be seen where there are red modifications and a lack of grey nucleosomes.

![alt text](https://github.com/cafelton/modbamvisualization/blob/main/nuc-threshold.png?raw=true)


