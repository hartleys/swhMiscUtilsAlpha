#!/usr/bin/env python
from __future__ import print_function

import sys,gzip,errno,time

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

eprint("STARTING");

args = sys.argv;

#+SHORTDESC="Truncates overlong cells so they don't break excel"
#+ HELPTEXT="..."$'\n'
#++HELPTEXT=" "$'\n'
#+ SYNTAXTEXT="truncateCells.py [options] infile.txt > outfile.txt"$'\n'
#+ PARAMS="--help: displays syntax and help info."$'\n'
#++PARAMS="--man: synonym for --help."$'\n'
#+ VERSION="0.0.5"


if "--help" in args or "--man" in args:
   import os;
   scriptAbsPath=os.path.abspath(__file__);
   shareMiscHome=os.path.dirname(os.path.abspath( os.path.dirname( __file__ )));
   os.system(shareMiscHome+"/internal.manUtil/manForScript "+scriptAbsPath);
   exit(1);

delim = "\t";
cellLim = 30000;
subDelim = -1;
overflowText = "CELL_LEN_OVERFLOW"
includeLen = False;

scriptString = args.pop(0);

while len(args) > 0:
   if "--delim" in args:
      idx = args.index("--delim");
      delim = args.pop(idx+1);
      args.pop(idx);
   elif "--cellLim" in args:
      idx = args.index("--cellLim");
      cellLim = int(args.pop(idx+1));
      args.pop(idx);
   elif "--subDelim" in args:
      idx = args.index("--subDelim");
      subDelim = args.pop(idx+1);
      args.pop(idx);
   elif "--overflowText" in args:
      idx = args.index("--overflowText");
      overflowText = args.pop(idx+1);
      args.pop(idx);
   elif "--includeLen" in args:
      idx = args.index("--includeLen");
      args.pop(idx);
      includeLen = True;
   else:
      break;

eprint("delim="+delim);
eprint("cellLim="+str(cellLim));
eprint("subDelim="+str(subDelim));
eprint("overflowText="+overflowText);
eprint("includeLen="+str(includeLen));


if len(args) == 0:
   inf = "-";
elif len(args) > 1:
   eprint("Ignoring unknown args: "+" ".join(args[0:-1]));
   inf = args.pop(-1);
else:
   inf = args.pop(-1);

if inf == "-":
   indata = sys.stdin;
elif inf[-3:] == ".gz":
   indata = gzip.open(inf,'r');
else:
   indata = open(inf,'r');


out = sys.stdout;


lnct = 0;
overflowLen = 0;
for line in indata:
   lnct = lnct + 1;
   if lnct % 1000 == 0:
      eprint(".", end='');
      if lnct % 5000 == 0:
         eprint(" ",end='');
      if lnct % 10000 == 0:
         eprint("["+str(lnct)+" lines complete] ["+time.strftime("%c")+"]")
   cells = line[0:-1].split(delim);
   firstCell = True;
   for cell in cells:
      if firstCell:
         texStart = "";
         firstCell = False;
      else:
         texStart = delim;
      if(subDelim == -1):
         if(len(cell) > cellLim):
            overflowLen = overflowLen + 1;
            if(includeLen):
               out.write(texStart+overflowText+"["+str(len(cell))+"]");
            else:
               out.write(texStart+overflowText);
         else:
            out.write(texStart+cell);
      else:
         subCells = cell.split(subDelim)
         #if(lnct < 20):
         #   eprint("subCells["+str(len(subCells))+"]="+"::".join(subCells));
         if(len(subCells) > cellLim):
            overflowLen = overflowLen + 1;
            if(includeLen):
               out.write(texStart+overflowText+"["+str(len(subCells))+"]");
            else:
               out.write(texStart+overflowText);
         else:
            out.write(texStart+cell);
   out.write("\n");


eprint("OverFlow Cells: "+str(overflowLen));
eprint("Total num Lines: " + str(lnct));