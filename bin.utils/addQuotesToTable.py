#!/usr/bin/env python

#TODO-WRITE-HELPDOC

from __future__ import print_function

#+SHORTDESC="Adds quotation marks to every entry of a table."
#+ HELPTEXT="Adds quotation marks to every entry of a table."$'\n'
#++HELPTEXT=""$'\n'
#+ SYNTAXTEXT="...."$'\n'
#++SYNTAXTEXT="...."$'\n'
#++SYNTAXTEXT="...."$'\n'
#+ PARAMS="--help or --man: displays syntax and help info."
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

#if "--help" in args or "--man" in args:
#   

outfile = args.pop(-1);
infile = args.pop(-1);


eprint("infile       = \""+infile+"\"");
eprint("outfile      = \""+outfile+"\"");

if "--tableDelimiter" in args:
   idx = args.index("--tableDelimiter");
   tableDelimiter = args.pop(idx+1);
   args.pop(idx);
else:
   tableDelimiter = "\t";

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

for line in indata:
   cells = line[:-1].split(tableDelimiter);
   outdata.write(tableDelimiter.join(["\""+s+"\"" for s in cells]) + "\n");

outdata.close();
