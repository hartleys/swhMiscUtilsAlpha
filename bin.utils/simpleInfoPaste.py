#!/usr/bin/env python
from __future__ import print_function

import sys,gzip,errno,time,random,collections

#+SHORTDESC="Simple paste-like merge of INFO fields in VCF"
#+ HELPTEXT="simpleInfoPaste.py utility. Warning: the two VCF files must"$'\n'
#++HELPTEXT="  be one to one and in the same exact order, or this tool will fail."$'\n'
#+ SYNTAXTEXT="simpleInfoPaste.py file1.vcf.gz file2.vcf.gz > outfile.vcf"$'\n'
#+ PARAMS="--help: displays syntax and help info."$'\n'
#++PARAMS="--man: displays syntax and help info."$'\n'
#+ VERSION="0.0.5"

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



if "--help" in args or "--man" in args:
   import os;
   scriptAbsPath=os.path.abspath(__file__);
   shareMiscHome=os.path.dirname(os.path.abspath( os.path.dirname( __file__ )));
   os.system(shareMiscHome+"/internal.manUtil/manForScript "+scriptAbsPath);
   exit(1);




INFILE1=args[1];
INFILE2=args[2];

out = sys.stdout;

inA=openSmart(INFILE1,"r");
inB=openSmart(INFILE2,"r");

lineA = inA.readline();
lineB = inB.readline();

infoSet=set();

while lineA.startswith("##"):
   out.write(lineA);
   if lineA.startswith("##INFO="):
      infokey = lineA.split(",")[0][11:]
      eprint("INFO=\""+infokey+"\"");
      infoSet.add( infokey )
   lineA = inA.readline();

while lineB.startswith("##"):
   if lineB.startswith("##INFO="):
      infokey = lineB.split(",")[0][11:]
      if infokey not in infoSet:
         out.write(lineB);
      else:
         eprint("WARNING: EXISTING TAG FOUND:"+infokey);
   lineB = inB.readline();

out.write(lineA);
lineA = inA.readline();
lineB = inB.readline();

lnct=0;

while lineA != "" and lineB != "":
   cellsA = lineA[:-1].split("\t");
   cellsB = lineB[:-1].split("\t");
   lnct = lnct + 1;
   if lnct % 1000 == 0:
         eprint(".", end='');
         if lnct % 5000 == 0:
            eprint(" ",end='');
         if lnct % 10000 == 0:
            eprint("["+str(lnct)+" lines complete] ["+time.strftime("%c")+"] ["+cellsA[0]+":"+cellsA[1]+" "+cellsA[3]+">"+cellsA[4]+"]")
   infoDelim=";";
   if cellsA[7][-1] == ";":
      infoDelim="";
   altAlleA = cellsA[4].split(",")[0];
   altAlleB = cellsB[4].split(",")[0];
   if cellsA[0] != cellsB[0] or cellsA[1] != cellsB[1] or cellsA[3] != cellsB[3] or altAlleA != altAlleB:
      eprint("ERROR ERROR ERROR: VCFs are not matched one-to-one!");
   cellsA[7] = cellsA[7] + infoDelim + cellsB[7];
   out.write("\t".join(cellsA)+"\n")
   lineA = inA.readline();
   lineB = inB.readline();


if lineA != "":
   eprint("ERROR ERROR ERROR: VCF files don't have the same length! file 1 is longer!");

if lineB != "":
   eprint("ERROR ERROR ERROR: VCF files don't have the same length! file 2 is longer!");




out.close();





