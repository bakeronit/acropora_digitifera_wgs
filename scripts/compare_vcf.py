#!/usr/bin/env python3

import argparse,re
import sys
import vcf
import collections
import pdb

parser = argparse.ArgumentParser(description='Compare pruned and imputed genotypes')
parser.add_argument('impvcf',metavar='imputed',type=argparse.FileType('r'),help='Imputed VCF File',default=sys.stdin)
parser.add_argument('pruned',metavar='pruned',type=argparse.FileType('r'),help='Pruned genotypes in text format',default=sys.stdin)
parser.add_argument('-o','--output',metavar='FILE',nargs='?',type=argparse.FileType('w'),help='Comparison output',default=sys.stdout)
parser.add_argument('-p',type=int,help='Num of sample to prune',default=1)

args  = parser.parse_args()
n_prune = args.p
vcf_reader = vcf.Reader(args.impvcf)

pruned_fields = args.pruned.readline().strip().split('\t')

#pdb.set_trace()

p_gt_template = collections.namedtuple('prunedgenotype',pruned_fields)

next_pruned = p_gt_template(*args.pruned.readline().strip().split('\t'))

args.output.write("record\tIMP_GT\tGT\tREF\tALT\n")

record_num = 1
count = 0
total = 1
heter = 0
t_heter = 0
for record in vcf_reader:
    for i in range(n_prune):
        if record_num == int(next_pruned.record):
            if next_pruned.GT in ['1/0','0/1','0|1','1|0']:
                t_heter += 1
            if collections.Counter(re.split('\/|\|',next_pruned.GT)) == collections.Counter(record.samples[int(next_pruned.sample_index)].data.GT.split("|")):
                count += 1
                if next_pruned.GT in ['1/0','0/1','0|1','1|0']:
                    heter += 1
            args.output.write(next_pruned.record+"\t"+record.samples[int(next_pruned.sample_index)].data.GT+"\t"+next_pruned.GT+"\t"+str(record.REF)+"\t"+str(record.ALT)+"\n")
            line = args.pruned.readline().strip()
            if line == "":
                break
            next_pruned = p_gt_template(*line.split('\t'))
            total += 1
    record_num+=1

print("The consistent rate under %d pruned is %.4f, at heterozygous site is %.4f"%(n_prune,count/total,heter/t_heter))

