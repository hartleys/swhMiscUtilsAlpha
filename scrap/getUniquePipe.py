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

if "--tableDelimiter" in args:
   idx = args.index("--tableDelimiter");
   tableDelimiter = args.pop(idx+1);
   args.pop(idx);
else:
   tableDelimiter = "\t";

if "--colNum" in args:
   idx = args.index("--colNum");
   colNum = int(args.pop(idx+1));
   args.pop(idx);
else:
   colNum = -1;

if infile == "-":
   indata = sys.stdin;
elif infile[-3:] == ".gz":
   indata = gzip.open(infile,'r');
else:
   indata = open(infile,'r');

outdata = sys.stdout;

outSet = set();

if colNum == -1:
   for line in indata:
      cell = line[:-1];
      if not cell in outSet:
         outdata.write(cell+"\n");
         outSet = outSet | set([cell]);
else:
   for line in indata:
      cell = line[:-1].split(tableDelimiter)[colNum];
      if not cell in outSet:
         outdata.write(cell+"\n");
         outSet = outSet | set([cell]);



