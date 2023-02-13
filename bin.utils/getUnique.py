#!/usr/bin/env python

from __future__ import print_function

#+SHORTDESC="Extracts/counts unique elements from a list"
#+ HELPTEXT="getUnique.py: Extracts and/or counts unique elements from a text list."$'\n'
#++HELPTEXT=" "$'\n'
#+ SYNTAXTEXT="getUnique.py [options] mylist.txt > myUniqueList.txt"$'\n'
#++SYNTAXTEXT="cat mylist.txt | getUnique.py [options] - > myUniqueList.txt"$'\n'
#+ PARAMS="--help: displays syntax and help info."$'\n'
#++PARAMS="--man: synonym for --help."$'\n'
#++PARAMS="--tableDelimiter <delim>: The field delimiter."$'\n'
#++PARAMS="--colNum <colnum>: Which column to take to fold uniques. By default this utility takes the entire line. If this option is specified, only the given column will be collapsed. All other columns will be discarded."$'\n'
#++PARAMS="--countReps: output an additional column with the number of times each unique element appears."$'\n'
#++PARAMS="--countCol <colnum>: The specified column should be composed of numeric values. If this param is used, output an additional column with the sum of this column for each unique value."$'\n'
#+ VERSION="0.0.5"


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



if len(args) == 0 or "--help" in args or "--man" in args:
   import os;
   scriptAbsPath=os.path.abspath(__file__);
   shareMiscHome=os.path.dirname(os.path.abspath( os.path.dirname( __file__ )));
   os.system(shareMiscHome+"/internal.manUtil/manForScript "+scriptAbsPath);
   exit(1);

#args = ["--colNum","0","../scrap/diseases.txt"]

infile = args.pop(-1);

if "--tableDelimiter" in args:
   idx = args.index("--tableDelimiter");
   tableDelimiter = args.pop(idx+1);
   args.pop(idx);
else:
   tableDelimiter = "\t";

if "--countReps" in args:
   countReps = True;
   args.pop(args.index("--countReps"));
else:
   countReps = False;

isSorted = False;
if "--isSorted" in args:
   isSorted = True;
   args.pop(args.index("--isSorted"));


if "--countCol" in args:
   idx = args.index("--countCol");
   countCol = int(args.pop(idx+1));
   args.pop(idx);
else:
   countCol = -1;

collapseCol=""
if "--collapseCol" in args:
   idx = args.index("--collapseCol");
   collapseCol = args.pop(idx+1);
   args.pop(idx);


#eprint("countCol = "+str(countCol));

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
outCountMap  = collections.defaultdict(lambda: 0);
outCountMap2 = collections.defaultdict(lambda: 0);

outCollapseCols = collections.defaultdict(lambda: frozenset());

eprint("starting...")

if colNum == -1 and countCol == -1 and (not countReps) and collapseCol == "":
  eprint("simple uniquify");
  for line in indata:
     cell = line[:-1];
     if cell not in outSet:
        outdata.write(cell+"\n");
     outSet.add(cell);
elif isSorted:
  # NOT FULLY TESTED!
  eprint("fast sorted uniquify");
  lnct = 0;
  prevct = 1;
  prevsum = 0;
  prevcell = "";
  for line in indata:
        if colNum == -1:
           cell = line[:-1];
        else:
           cell = line[:-1].split(tableDelimiter)[colNum];
        if cell != prevcell and lnct != 0:
           countStr1 = "" if not countReps else tableDelimiter+str(prevct);
           countStr2 = "" if countCol == -1 else tableDelimiter+str(prevsum);
           outdata.write( prevcell+countStr1+countStr2+"\n" )
           prevcell = cell;
           prevct = 0;
           prevsum = 0;
        elif cell != prevcell:
           prevcell = cell;
           prevct = 0;
           prevsum = 0;
        if countCol >= 0:
           prevsum = prevsum + int(line[:-1].split(tableDelimiter)[countCol]);
        prevct = prevct + 1;
        lnct = lnct + 1;
  countStr1 = "" if not countReps else tableDelimiter+str(prevct);
  countStr2 = "" if countCol == -1 else tableDelimiter+str(prevsum);
  outdata.write( prevcell+countStr1+countStr2+"\n" )
else:
  eprint("complex uniquify");
  for line in indata:
        if colNum == -1:
           cell = line[:-1];
        else:
           cell = line[:-1].split(tableDelimiter)[colNum];
        outSet = outSet | set([cell]);
        outCountMap[cell] = outCountMap[cell] +1;
        if countCol >= 0:
           outCountMap2[cell] = outCountMap2[cell] + int(line[:-1].split(tableDelimiter)[countCol]);

  for cell in sorted(outSet):
     countStr1 = "" if not countReps else tableDelimiter+str(outCountMap[cell]);
     countStr2 = "" if countCol == -1 else tableDelimiter+str(outCountMap2[cell]);
     outdata.write(cell + countStr1 + countStr2 + "\n");

