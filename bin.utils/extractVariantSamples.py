#!/usr/bin/env python

from __future__ import print_function

import sys,gzip,errno,collections

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

def is_integer(s):
    try:
       int(s);
       return True;
    except ValueError:
       pass
    return False;


args = sys.argv;

#args = ["--colNum","0","../scrap/diseases.txt"]

infile = args.pop(-1);


if infile == "-":
   indata = sys.stdin;
elif infile[-3:] == ".gz":
   indata = gzip.open(infile,'r');
else:
   indata = open(infile,'r');

outdata = sys.stdout;

#outCountMap2 = collections.defaultdict(lambda: 0);

sampleDict = dict();

for line in indata:
   if line[0] == "#":
      if line[1] != "#":
         #read in sample ids:
         i = 0;
         cells = line[:-1].split("\t");
         for samp in cells[9:]:
            sampleDict[i] = samp;
            i = i + 1;
         outdata.write("\t".join(cells[0:8]) +"\tVariantSamples...\n")
   else:
      #read in vcf lines:
      cells = line[:-1].split("\t");
      outCells = [];
      i = 0;
      for geno in cells[9:]:
         if (geno[0] != "0" and geno[0] != ".") or (geno[2] != "0" and geno[2] != "."):
            outCells = outCells + [sampleDict[i] + ":" + geno];
         i = i + 1;
      outdata.write("\t".join(cells[0:8])+"\t"+"\t".join(outCells)+"\n");

outdata.close();
#
#
#line = open("testtest.txt",'r').readline();
#for geno in cells[9:]:
#   if (geno[0] != "0" and geno[0] != ".") or (geno[2] != "0" and geno[2] != "."):
#      outCells = outCells + [geno];
#   i = i + 1;
#
#