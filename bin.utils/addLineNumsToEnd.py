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
#+ HELPTEXT="vcf2table.py utility. Takes the info data from a VCF and makes a table. For testing purposes."$'\n'
#++HELPTEXT=" "$'\n'
#+ SYNTAXTEXT="vcf2table.py [options] myVcf.vcf.gz"$'\n'
#+ PARAMS="--help: displays syntax and help info."$'\n'
#++PARAMS="--man: synonym for --help."$'\n'
#++PARAMS="--infoColumnList <tag1,tag2,...>: a comma-delimited list of INFO columns you want extracted, in order."$'\n'
#++PARAMS="--infoColumnReordering <tag1,tag2,...>: a comma-delimited list of INFO columns you want to be placed to the left of the other info columns, in order."$'\n'
#++PARAMS="--keepGT: Flag, keep genotype columns."$'\n'
#++PARAMS="--filterINFO [LIST]: ..."$'\n'
#++PARAMS="--filterGT [LIST]: ..."$'\n'
#++PARAMS="--filterFilter: If set, drop variants that fail the FILTER column."$'\n'
#++PARAMS="--tableGT: ..."$'\n'
#++PARAMS="--flatTableGT: ..."$'\n'
#++PARAMS="--gtLimit: ..."$'\n'
#++PARAMS="--gtColumnList: ..."$'\n'
#++PARAMS="--gtTag: ..."$'\n'
#++PARAMS="--gtDelim: ..."$'\n'
#++PARAMS="--emitOnlyNonref: ..."$'\n'
#++PARAMS="--emitSampleIdsOnly: ..."$'\n'
#++PARAMS="--infoIsNumericList [LIST]: Comma-delimited list of INFO tags that should be formatted as numerics, as opposed to quoted."$'\n'
#++PARAMS="--headerInfoFile <filename>: If this parameter is used, then a file will be created containing the INFO column names and descriptions from the header, in the same order as the columns in the output table."$'\n'
#++PARAMS="--extendedHeaderInfoFile <filename>: If this parameter is used, then an extended header file will be created with extra information."$'\n'
#++PARAMS="--sampleDecoder <filename>: Used with the tableGT option. This provides a sample-level table with phenotype information which will be merged into the genotype table. The first column must be sample.ID"$'\n'
#++PARAMS="--closeQuote: ..."$'\n'
#++PARAMS="--openQuote: ..."$'\n'
#++PARAMS="--missingVal [value]: The string to be used to indicate missing values"$'\n'
#++PARAMS="--noQuote: Do not add excel-formatted quotes."$'\n'
#++PARAMS="--replaceUnderscores: In all column names, replace all underscores with spaces."$'\n'
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


eprint("INDATA="+args[-1]);

if args[-1] == "-":
   args.pop(-1)
   inf = sys.stdin;
   isStream = True;
else:
   inf = args.pop(-1);
   isStream = False;


if not isStream:
   if inf[-3:] == ".gz":
      inf = gzip.open(inf,'r');
   else:
      inf = open(inf,'r');

MAXSPAN=100

out = sys.stdout;



lnct = 0;
for line in inf:
    lnct = lnct + 1;
    cells = line[:-1].split("\t");
    if lnct % 10000 == 0:
       eprint(".", end='');
       if lnct % 50000 == 0:
          eprint(" ",end='');
       if lnct % 100000 == 0:
          eprint("["+str(lnct)+" lines complete] ["+time.strftime("%c")+"] ["+line[:-1]+"]")
    out.write(line[:-1]+"\t"+str(lnct)+"\n");
out.close();


