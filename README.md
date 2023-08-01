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

