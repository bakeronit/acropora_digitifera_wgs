#!/usr/bin/env python

"""
MIT License

Copyright (c) 2020 Michael Alonge <malonge11@gmail.com>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""

import os
import sys
import argparse
from collections import defaultdict

from AGPFile import AGPFile


def sup_update(in_file, agp_file, args):
    # Make a dictionary associating each original sequence with the destination sequence
    trans = {}
    strands = {}
    seq_lens = {}
    agp = AGPFile(agp_file, mode="r")
    for agp_line in agp.iterate_lines():
        if not agp_line.is_gap:
            start, end = agp_line.obj_beg - 1, agp_line.obj_end
            trans[agp_line.comp] = (start, end, agp_line.obj)
            strands[agp_line.comp] = agp_line.orientation
            seq_lens[agp_line.comp] = end - start

    # Iterate through the gff intervals and update them according to trans
    with open(in_file, "r") as f:
        for line in f:
            line = line.rstrip()
#            print(line)
            if line.startswith("#"):
                print(line)  # Print this comment line
            else:
                fields = line.split("\t")
                if args.file_type=='gff':
                    h, s, e, st = fields[0], int(fields[3]), int(fields[4]), fields[6]
                    s -= 1  # Keep everything zero-indexed
                elif args.file_type=='bed':
                    h, s, e, st = fields[0], int(fields[1]), int(fields[2]), fields[4]
                elif args.file_type=='pos':
                    h, s = fields[0], int(fields[1])
                    e=s
                    st="+"
                else:
                    raise Error("Invalid args.file_type")

                if h not in trans:
                    sys.stderr.write(line)
                    continue
#                    raise ValueError("Inconsistent input files.")
                l = seq_lens[h]
                
                if e > l:
                    sys.stderr.write("Excluded out of bounds feature "+line+"\n")
                    continue

                # Check if the original sequence has been reverse complemented
                if strands[h] == "-":

                    s, e = l-e, l-s
                    if st == "+":
                        st = "-"
                    else:
                        st = "+"


                new_s = trans[h][0] + s
                fields[0] = trans[h][2]
                if ( args.file_type=='gff'):
                    new_e = trans[h][0] + e
                    fields[3] = str(new_s + 1)  # back to one-based indexing for gff format
                    fields[4] = str(new_e)
                    fields[6] = st
                elif ( args.file_type=='bed'):
                    new_e = trans[h][0] + e
                    fields[1] = str(new_s)  
                    fields[2] = str(new_e)
                    fields[4] = st
                elif ( args.file_type=='pos'):
                    fields[1] = str(new_s)

                if args.keep:
                    print(line+ "\t"+"\t".join(fields))
                else:
                    print("\t".join(fields))



def main():
    parser = argparse.ArgumentParser(description="Update gff intervals given a RagTag AGP file", usage="ragtag.py updategff [-c] <genes.gff> <ragtag.agp>")

    parser.add_argument("coordinate_file", metavar="FILE", default="", type=str, help="gff/bed/pos file")
    
    parser.add_argument('--file-type',
    dest='file_type',
    choices=['gff','bed','pos'],
    default='pos',
    help='Filetype for coordinates to be translated (Default: pos)')

    parser.add_argument('--keep',
    action='store_true',    
    dest='keep',
    help='Print the original (old) coordinates as well as the new ones')

    parser.add_argument("agp", nargs='?', default="", metavar="<ragtag.*.agp>", type=str, help="agp file")

    args = parser.parse_args()

    if not args.coordinate_file or not args.agp:
        parser.print_help()
        sys.exit()

    coord_file = os.path.abspath(args.coordinate_file)
    agp_file = os.path.abspath(args.agp)


    sup_update(coord_file, agp_file, args)



#    log("Goodbye")


if __name__ == "__main__":
    main()
