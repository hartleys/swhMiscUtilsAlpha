#!/usr/bin/env python

from __future__ import print_function

#TODO: WRITE HELPDOC?

import sys,gzip,errno,collections

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

args = sys.argv;

infileA = args[-6];
fromColA = int(args[-5]);
toColA = int(args[-4]);
infileB = args[-3];
fromCol = int(args[-2]);
outfile = args[-1];

#args = ["entrezToEnsembl.txt","0","1","ensemblToKnown.txt","0","1","entrezToEnsemblToKnown.txt"]


eprint("decoderFile  = \""+infileA+"\"");
eprint("decFromCol = \""+str(fromColA)+"\"");
eprint("decToCol   = \""+str(toColA)+"\"");
eprint("infile  = \""+infileB+"\"");
eprint("fromCol = \""+str(fromCol)+"\"");
eprint("outfile  = \""+outfile+"\"");

if infileA == "-":
   indataA = sys.stdin;
elif infileA[-3:] == ".gz":
   indataA = gzip.open(infileA,'r');
else:
   indataA = open(infileA,'r');

if infileB == "-":
   indataB = sys.stdin;
elif infileB[-3:] == ".gz":
   indataB = gzip.open(infileB,'r');
else:
   indataB = open(infileB,'r');

if outfile == "-":
   outdata = sys.stdout;
elif outfile[-3:] == ".gz":
   outdata = gzip.open(outfile,'w');
else:
   outdata = open(outfile,'w');

fromDictA = collections.defaultdict(lambda: set());
for line in indataA:
   cells = line[:-1].split('\t');
   if cells[toColA] != "" and cells[fromColA] != "":
      for fr in cells[fromColA].split(","):
         fromDictA[fr] = fromDictA[fr] | set(cells[toColA].split(","));

countTableTemp = collections.defaultdict(lambda: 0);
for key in sorted(fromDictA.keys()):
   countTableTemp[len(fromDictA[key])] = countTableTemp[len(fromDictA[key])] + 1;

eprint("Dict A:");

for ct in sorted(countTableTemp.keys()):
   eprint(str(ct)+"\t"+str(countTableTemp[ct]));

for line in indataB:
   cells = line[:-1].split('\t');
   outColString = "";
   if cells[fromCol] != "":
      currSet = set();
      for fr in cells[fromCol].split(","):
         currSet = currSet | fromDictA[fr];
      if(len(currSet) > 0):
         outColString = ",".join(sorted(currSet));
   outdata.write(line[:-1] + "\t"+outColString + "\n");

outdata.close();



