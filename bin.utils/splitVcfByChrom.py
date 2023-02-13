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


#+SHORTDESC="Splits a VCF into parts by chromosome."
#+ HELPTEXT="splitVcfByChrom.py: Splits a VCF into parts by chromosome."$'\n'
#++HELPTEXT=" "$'\n'
#+ SYNTAXTEXT="cat myVCF.vcf | splitVcfByChrom.py [options] output.prefix"$'\n'
#+ PARAMS="--help: displays syntax and help info."$'\n'
#++PARAMS="--man: synonym for --help."$'\n'
#+ VERSION="0.0.5"


if "--help" in args or "--man" in args:
   import os;
   scriptAbsPath=os.path.abspath(__file__);
   shareMiscHome=os.path.dirname(os.path.abspath( os.path.dirname( __file__ )));
   os.system(shareMiscHome+"/internal.manUtil/manForScript "+scriptAbsPath);
   exit(1);



#args = ["--colNum","0","../scrap/diseases.txt"]


indata = sys.stdin;
outPrefix = args.pop(-1);

headerLines = [];
currChrom = "";
#out = "";

for line in indata:
   if line[0] == "#":
      headerLines = headerLines + [line];
   else:
      cells = line.split("\t",1);
      if cells[0] != currChrom:
         eprint("starting chrom: " + cells[0]);
         out = open(outPrefix + "." + cells[0] + ".vcf",'w');
         currChrom = cells[0];
         for hl in headerLines:
            out.write(hl);
      out.write(line);





