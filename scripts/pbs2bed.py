#!/usr/bin/env python3

import argparse
import sys
import math

parser = argparse.ArgumentParser(description='Identify high PBS regions and print in bed format')

parser.add_argument('inpbs',metavar='FILE',nargs='?',type=argparse.FileType('r'),help='Input PBS File',default=sys.stdin)
parser.add_argument('-t','--pbs-threshold',type=float,help='Keep regions with PBS above this value',default=20)


args  = parser.parse_args()

start_pos = None
end_pos = None
current_scaff = None
scores = []

for line in args.inpbs:

	line_values = line.strip().split()
	scaff = str(line_values[0])
	pbs = float(line_values[2])
	if not math.isfinite(pbs): pbs = 5

	if current_scaff==None:
		current_scaff=scaff

#	sys.stdout.write("l: "+line.strip()+" lv:\t"+line_values[2]+":\t"+str(pbs)+"\n")
	if (pbs > args.pbs_threshold) and (scaff==current_scaff):
		# Either start a new record or continue an old one on the same scaffold
		if start_pos == None:
			start_pos = float(line_values[1])

		end_pos = float(line_values[1])
		scores.append(pbs)

	else:

		# Write a valid record if we have one
		if ( start_pos != None ):
			# This means we have a valid record to write
			sys.stdout.write(scaff+"\t"+str(math.floor(start_pos))+"\t"+
										str(math.floor(end_pos))+"\t"+str(round(max(scores),1))+"\t"+str(round(sum(scores),1))+"\n")
		if (scaff!=current_scaff and pbs > args.pbs_threshold):
			start_pos = float(line_values[1])
			end_pos = start_pos
			scores = [pbs]
			current_scaff=scaff
		else:
			start_pos = None
			scores=[]
