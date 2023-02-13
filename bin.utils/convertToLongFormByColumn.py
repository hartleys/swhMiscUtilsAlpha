#!/usr/bin/env python
from __future__ import print_function

import sys,gzip,errno,time,random,collections

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)
#TODO: WRITE HELPDOC?

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

def dexcel(s):
    #If a string has single or double quotes around it, remove them.
    #Make sure the pair of quotes match.
    #If a matching pair of quotes is not found, return the string unchanged.
    ss = s;
    if s.startswith("="):
        ss = s[1:];
    if (len(ss) >= 2) and (ss[0] == ss[-1]) and (ss.startswith(("'", '"'))):
        return ss[1:-1]
    return ss

eprint("STARTING");

args = sys.argv;

#+SHORTDESC="Converts VCF to tab-delimited table"
#+ HELPTEXT="vcf2table.py utility. Takes the info data from a VCF and makes a table. For testing purposes."$'\n'
#++HELPTEXT=" "$'\n'
#+ SYNTAXTEXT="vcf2table.py [options] myVcf.vcf.gz"$'\n'
#+ PARAMS="--help: displays syntax and help info."$'\n'
#++PARAMS="--man: synonym for --help."$'\n'
#+ VERSION="0.0.5"


if len(args) == 0 or "--help" in args or "--man" in args:
   import os;
   scriptAbsPath=os.path.abspath(__file__);
   shareMiscHome=os.path.dirname(os.path.abspath( os.path.dirname( __file__ )));
   os.system(shareMiscHome+"/internal.manUtil/manForScript "+scriptAbsPath);
   exit(1);



eprint("INDATA="+args[-1]);




NULLCHAR="-"
if "--missingVal" in args:
   idx = args.index("--missingVal");
   args.pop(idx);
   NULLCHAR=args.pop(idx)

TOOMANYSTRING="-"
if "--tooManyString" in args:
   idx = args.index("--tooManyString");
   args.pop(idx);
   TOOMANYSTRING=args.pop(idx)


maxct=30
if "--maxCt" in args:
   idx = args.index("--maxCt");
   args.pop(idx);
   maxct=int(args.pop(idx))


hasNoSplitNames=True
splitNames=[];
if "--splitNames" in args:
   idx = args.index("--splitNames");
   args.pop(idx);
   splitNames=args.pop(idx).split(",");
   hasNoSplitNames=False;

masterSplit="";
if "--master" in args:
   idx = args.index("--master");
   args.pop(idx);
   masterSplit=args.pop(idx)

splitColumnNameString=args.pop(1);
splitColumnNames=splitColumnNameString.split(",");
if hasNoSplitNames:
   splitNames = splitColumnNames;

for ii in range(0,len(splitColumnNames)):
   scn = splitColumnNames[ii];
   sn = splitNames[ii];
   eprint("splitting on: "+scn + " / "+sn);


if args[-1] == "-" or len(args) == 1:
   args.pop(-1)
   inf = sys.stdin;
   isStream = True;
   eprint("reading from stdin");
else:
   inf = args.pop(-1);
   isStream = False;
   if inf[-3:] == ".gz":
      eprint("reading gzip file: "+inf);
      inf = gzip.open(inf,'r');
   else:
      eprint("reading text file: "+inf);
      inf = open(inf,'r');

out = sys.stdout;


eprint("Starting Iteration ["+time.strftime("%c")+"]")

lnct=0
splitIdx = [-1 for i in splitColumnNames]
masterIdx = -1;

isTooManyCt = 0;

try:
   for line in inf:
      lnct = lnct + 1;
      if lnct % 1000 == 0 :
         eprint(".", end='');
         if lnct % 5000 == 0:
            eprint(" ",end='');
         if lnct % 10000 == 0:
            cells = line.split("\t")
            eprint("["+str(lnct)+" lines complete] ["+time.strftime("%c")+"]")
      if lnct == 1:
         eprint("parsing title line...");
         cells = line[:-1].split("\t");
         for i in range(0,len(cells)):
            c = cells[i];
            eprint("c = ",c);
            if c in splitColumnNames:
               splitIdx[splitColumnNames.index(c)] = i;
               eprint("splitColumnNames["+str(splitColumnNames.index(c))+"]="+str(splitIdx[splitColumnNames.index(c)]));
            if c == masterSplit and masterSplit != "":
               masterIdx = i;
         out.write(line[:-1]+"\tTYPE\tSAMP\n");
      else:
         cells = [dexcel(x) for x in (line[:-1].split("\t"))];
         outputStrings = [];
         outputTypes = [];
         foundTooMany = False;
         for ii in range(0,len(splitIdx)):
            sidx = splitIdx[ii];
            if cells[sidx] != "." and cells[sidx] != "-" and cells[sidx] != "TOO_MANY_TO_PRINT":
               ccells = cells[sidx].split(",");
               for cc in ccells:
                  #out.write(line[:-1]+"\t"+splitColumnNames[ii]+"\t"+cc+"\n");
                  outputStrings = outputStrings + [cc]
                  outputTypes = outputTypes + [splitNames[ii]]
            elif cells[sidx] == "TOO_MANY_TO_PRINT":
               foundTooMany = True;
            #   out.write(line[:-1]+"\t"+splitColumnNames[ii]+"\tTOO_MANY_TO_PRINT\n");
         if (not foundTooMany) and  len(outputStrings) > 0 and len(outputStrings) <= maxct:
            for ii in range(0,len(outputStrings)):
               out.write(line[:-1]+"\t"+outputTypes[ii]+"\t"+outputStrings[ii]+"\n");
         elif foundTooMany or len(outputStrings) > maxct:
            out.write(line[:-1]+"\t"+"TOOMANYTOPRINT"+"\t"+"TOOMANYTOPRINT"+"\n");
            isTooManyCt =  isTooManyCt +1;

except IOError as e:
   if e.errno == errno.EPIPE:
      exit(0);

eprint("processed ",str(lnct), " lines");
eprint("too many: ",str(isTooManyCt), "");


out.close();

