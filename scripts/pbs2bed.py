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
scores = []

for line in args.inpbs:

	line_values = line.split()
	scaff = str(line_values[0])
	pbs = float(line_values[2])
	if ( pbs > args.pbs_threshold ):
		if start_pos == None:
			start_pos = float(line_values[1])

		end_pos = float(line_values[1])
		scores.append(pbs)

	else:
		if ( start_pos != None ):
			sys.stdout.write(scaff+"\t"+str(math.floor(start_pos))+"\t"+
										str(math.floor(end_pos))+"\t"+str(round(max(scores),1))+"\n")
		start_pos = None
		scores=[]
