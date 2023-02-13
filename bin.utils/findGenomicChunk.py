#!/usr/bin/env python
from __future__ import print_function

import sys,gzip,errno,time,random,collections,re

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

def isNumber(s):
    try:
        float(s)
        return True
    except ValueError:
        pass
    return False

def intOrZero(s):
    try:
        return int(s)
    except ValueError:
        pass
    return 0

def floatOrZero(s):
    try:
        return float(s)
    except ValueError:
        pass
    return 0

def openSmart(s,t):
   if(s[-3:] == ".gz"):
      return gzip.open(s,t);
   else:
      return open(s,t);

def indexOf(lst,x):
    try:
       y = lst.index(x);
       return y;
    except ValueError:
       return -1;

def dequote(s):
    #If a string has single or double quotes around it, remove them.
    #Make sure the pair of quotes match.
    #If a matching pair of quotes is not found, return the string unchanged.
    if (s[0] == s[-1]) and s.startswith(("'", '"')):
        return s[1:-1]
    return s


args = sys.argv;

#+SHORTDESC="Converts VCF to tab-delimited table"
#+ HELPTEXT="findGenomicChunk.py utility."$'\n'
#++HELPTEXT=" "$'\n'
#+ SYNTAXTEXT="findGenomicChunk.py [options] /path/to/genomic/chunk/dir chr12 1242342"$'\n'
#++SYNTAXTEXT="cat genomicChunkList.txt | findGenomicChunk.py [options] /path/to/genomic/chunk/dir -"$'\n'
#++SYNTAXTEXT="   On CCAD: findGenomicChunk.py /mnt/nfs/gigantor/ifs/Shared/hartleys/resources/makeSplitIVs/split/p200.v01.win20k chr12 1242342"$'\n'
#++SYNTAXTEXT="   On CCAD: findGenomicChunk.py /data/hartleys/pub/resources/makeSplitIvs/split/p200.v01.win20k chr12 1242342"$'\n'
#+ PARAMS="--help: displays syntax and help info."$'\n'
#++PARAMS="--man: synonym for --help."$'\n'
#+ VERSION="0.0.5"


if len(args) == 0 or "--help" in args or "--man" in args:
   import os;
   scriptAbsPath=os.path.abspath(__file__);
   shareMiscHome=os.path.dirname(os.path.abspath( os.path.dirname( __file__ )));
   os.system(shareMiscHome+"/internal.manUtil/manForScript "+scriptAbsPath);
   exit(1);


if args[-1] == "-":
   ivDIR = args[-2];
   inf = sys.stdin;
   for chrpos in inf:
      chrom = chrpos[:-1].split("\t")[0];
      pos = int(chrpos[:-1].split("\t")[1]);
      hasAnyChunk = False
      for line in open(ivDIR+"/chunkInfo.txt",'r'):
         cells = line[:-1].split("\t");
         if cells[0] == chrom and pos >= int(cells[1]) and pos < int(cells[2]):
             sys.stdout.write(cells[3] + "\n");
             hasAnyChunk=True;
      if not hasAnyChunk:
         sys.stdout.write("NA\n");
else:
   ivDIR=args[-3];
   chrom=args[-2];
   pos = int(args[-1]);

   #eprint( "searching for: " + chrom+":"+str(pos));
   hasAnyChunk = False
   for line in open(ivDIR+"/chunkInfo.txt",'r'):
      cells = line[:-1].split("\t");
      #eprint("checking: "+line[:-1]);
      #eprint("checking: "+cells[0]+":"+cells[1]+"-"+cells[2]);
      if cells[0] == chrom and pos >= int(cells[1]) and pos < int(cells[2]):
          sys.stdout.write(cells[3] + "\n");
          hasAnyChunk=True;
   if not hasAnyChunk:
      sys.stdout.write("NA\n");


