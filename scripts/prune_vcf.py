#!/usr/bin/env python3

import argparse
import sys
import vcf
import random
import pdb
import collections

parser = argparse.ArgumentParser(description='Randomly remove genotypes from a vcf')
parser.add_argument('invcf',metavar='FILE',nargs='?',type=argparse.FileType('r'),help='Input VCF File',default=sys.stdin)
parser.add_argument('-o','--output',metavar='FILE',nargs='?',type=argparse.FileType('w'),help='Output VCF File',default=sys.stdout)
parser.add_argument('--save-filtered',metavar='FILE',type=argparse.FileType('w'),help='Save original genotypes to file',default=None)
parser.add_argument('--sites-file',metavar='FILE',type=argparse.FileType('r'),help='File listing genotypes to sample',default=None)
parser.add_argument('-n',type=int,help='Num site to prune',default=10000)
parser.add_argument('-l',type=int,help='Total num sites',required=True)
parser.add_argument('-p',type=int,help='Num of sample to prune',default=1)

args  = parser.parse_args()

vcf_reader = vcf.Reader(args.invcf)

vcf_writer = vcf.Writer(args.output,vcf_reader)

if args.save_filtered != None:
    args.save_filtered.write("record\tCHROM\tPOS\tREF\tALT\tsample_index\tGQ\tGT\n")
n_prune = args.p
n_sample = len(vcf_reader.samples)
next_site=None
filter_sites=None

if args.sites_file:
    fields = args.sites_file.readline().strip().split('\t')
    template = collections.namedtuple('prunedgenotype',fields)
    next_site = template(*args.sites_file.readline().strip().split('\t'))
else:
    filter_sites = set(random.sample(range(1, args.l), args.n))

def nullify_call(call):
    # Create blank call using original call format
    blank_call_data = [None] * len(call.data)
    blank_call_data[call.data._fields.index('GT')] = './.'

    call.data = call.data.__class__(*blank_call_data)
    return call

def should_sample(record_num,filter_sites,next_site):
    if filter_sites!=None:
        return (record_num in filter_sites)
    else:
        return (record_num == int(next_site.record))

record_num = 1
for record in vcf_reader:
    if should_sample(record_num,filter_sites,next_site):
        if filter_sites!=None:
            #filt_index=random.randint(0, n_sample - 1)
            filt_index_list = random.sample(range(0,n_sample),n_prune)
        else:
            filt_index=int(next_site.sample_index)
            next_site = template(*args.sites_file.readline().strip().split('\t'))
        
        for filt_index in filt_index_list:
            filt_call = record.samples[filt_index]
            original_gq = filt_call.data.GQ

            if args.save_filtered != None:
                args.save_filtered.write(str(record_num)+"\t"+record.CHROM+"\t"+str(record.POS)+"\t"+str(record.REF)+"\t"+str(record.ALT)+"\t"+str(filt_index)+"\t"+str(original_gq)+"\t"+ str(filt_call.data.GT)+"\n")
                args.save_filtered.flush()

            filt_call = nullify_call(filt_call)

            sys.stderr.write("Filtered genotype with GQ:"+str(original_gq)+" at sample "+str(filt_index)+" record "+str(record_num)+"\n")
            sys.stderr.flush()
    record_num+=1
    vcf_writer.write_record(record)
