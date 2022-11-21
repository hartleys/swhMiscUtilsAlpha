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


eprint("STARTING");

args = sys.argv;

#+SHORTDESC="Converts VCF to tab-delimited table"
#+ HELPTEXT="findGenomicChunk.py utility."$'\n'
#++HELPTEXT=" "$'\n'
#+ SYNTAXTEXT="findGenomicChunk.py [options] /path/to/genomic/chunk/dir chr12 1242342"$'\n'
#+ PARAMS="--help: displays syntax and help info."$'\n'
#++PARAMS="--man: synonym for --help."$'\n'
#+ VERSION="0.0.5"


if "--help" in args or "--man" in args:
   import inspect,commands,os
   SCRIPT_FILE_PATH=inspect.stack()[0][1];
   print(SCRIPT_FILE_PATH);
   SCRIPT_DIR_PATH=commands.getstatusoutput("dirname "+SCRIPT_FILE_PATH)[1]
   SCRIPT_FILENAME=commands.getstatusoutput("basename "+SCRIPT_FILE_PATH)[1]
   #print("SCRIPT_DIR_PATH="+SCRIPT_DIR_PATH);
   #print("SCRIPT_FILENAME="+SCRIPT_FILENAME);
   os.system(SCRIPT_DIR_PATH+"/../internal.manUtil/manForScript"+" "+SCRIPT_FILE_PATH)
   exit(0);


if args[-1] == "-":
   ivDIR = args[-2];
   inf = sys.stdin;
   for chrpos in inf:
      chrom = chrpos[-1].split("\t")[0];
      pos = int(chrpos[-1].split("\t")[1]);
      hasAnyChunk = False
      for line in open(ivDIR+"/chunkInfo.txt",'r'):
         cells = line[-1].split("\t");
         if cells[0] == chrom and pos >= int(cells[1]) and pos < int(cells[2]):
             sys.stdout.write(cells[3] + "\n");
             hasAnyChunk=True;
      if not hasAnyChunk:
         sys.stdout.write("NA\n");
else:
   ivDIR=args[-3];
   chrom=args[-2];
   pos = int(args[-1]);
   hasAnyChunk = False
   for line in open(ivDIR+"/chunkInfo.txt",'r'):
      cells = line[-1].split("\t");
      if cells[0] == chrom and pos >= int(cells[1]) and pos < int(cells[2]):
          sys.stdout.write(cells[3] + "\n");
          hasAnyChunk=True;
   if not hasAnyChunk:
      sys.stdout.write("NA\n");


