#!/usr/bin/env python

from __future__ import print_function

import sys,gzip,errno,collections

#print("???");

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

if "--tableDelimiter" in args:
   idx = args.index("--tableDelimiter");
   tableDelimiter = args.pop(idx+1);
   args.pop(idx);
else:
   tableDelimiter = "\t";

outfile = args.pop(-1);
dropCols = args.pop(-1).split(",");
infile = args.pop(-1);

#eprint("infile="+infile);
#eprint("dropCols="+"\t".join([str(d) for d in dropCols]));
#eprint("outfile="+outfile);

if infile == "-":
   indata = sys.stdin;
elif infile[-3:] == ".gz":
   indata = gzip.open(infile,'r');
else:
   indata = open(infile,'r');

if outfile == "-":
   outdata = sys.stdout;
elif outfile[-3:] == ".gz":
   outdata = gzip.open(outfile,'w');
else:
   outdata = open(outfile,'w');

if is_integer(dropCols[0]):
   dropCols = [int(d) for d in dropCols];
else:
   titleLine = indata.readline()[:-1];
   titleCells = titleLine.split(tableDelimiter);
   dropColNames = dropCols;
   dropCols = [i for i in range(len(titleCells)) if (titleCells[i] in dropColNames)];
   outdata.write(tableDelimiter.join([titleCells[i] for i in range(len(titleCells)) if not i in dropCols ]) + "\n");
   #eprint("filtering column names: " + ",".join(dropColNames));

#eprint("dropCols="+"\t".join([str(d) for d in dropCols]));

for line in indata:
   cells = line.split(tableDelimiter);
   outdata.write(tableDelimiter.join([cells[i] for i in range(len(cells)) if not i in dropCols ]));

#outdata.close();






