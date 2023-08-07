import sys


out = open('.'.join(sys.argv[1].split('.')[:-1]) + '-genes.bed', 'w')
for line in open(sys.argv[1]):
    if line[0] != '#':
        line = line.split('\t')
        if line[2] == 'gene':
            if 'Name=' in line[8]:
                genename = line[8].split('Name=')[1].split(';')[0]
            else: genename = line[8].split('gene_id=')[1].split(';')[0]
            out.write('\t'.join([line[0], line[3], line[4], genename, '0', line[6]]) + '\n')
out.close()



