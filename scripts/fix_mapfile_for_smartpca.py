#!/usr/bin/env python
import subprocess
from collections import defaultdict

reference = "reference.fa"
mapfile = 'Adigitifera_ldpruned.map'
faidx_cmd = f'faidx {reference} -i chromsizes'
faidx_out = subprocess.run(faidx_cmd, shell=True, stdout=subprocess.PIPE).stdout.splitlines()

cumlen = 0
cumlen_at_chrom = defaultdict(int)
for line in faidx_out:
    line = line.decode()
    chrom = line.split()[0]
    cumlen_at_chrom[chrom] = cumlen
    length = int(line.split()[1])
    cumlen += length

n = 1
with open(mapfile,'rt') as fh:
    for line in fh:
        line = line.strip()
        chrom = line.split()[0]
        pos = int(line.split()[-1])
        new_pos = cumlen_at_chrom[chrom] + pos
        print(f'1 rs{n}:{str(new_pos)} 0 {str(new_pos)}')
        n += 1 
