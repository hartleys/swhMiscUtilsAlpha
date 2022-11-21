#!/usr/bin/env python
from __future__ import print_function

import sys,gzip,errno,time,random,collections

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

#+SHORTDESC="Converts VCF to intervar input table"
#+ HELPTEXT="vcf2intervar.py utility."$'\n'
#++HELPTEXT=" "$'\n'
#+ SYNTAXTEXT="vcf2intervar.py [options] myVcf.vcf.gz > myIntervarTable.txt"$'\n'
#+ PARAMS="--help: displays syntax and help info."$'\n'
#++PARAMS="--man: synonym for --help."$'\n'
#++PARAMS="--skipLines n: skip n lines. IF your file has a title line, set this to 1. Note that VCF header lines will be skipped automatically."$'\n'
#++PARAMS="--keepMulti: Keep multiallelics. By default only the first allele will be kept. This is useful because varmyknife-multiallelic-split files will needlessly confuse intervar"$'\n'
#++PARAMS="--writeTitleLine: Write a title line. Not sure if intervar wants this or not."$'\n'
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


skipLines=0;
if "--skipLines" in args:
   idx = args.index("--skipLines");
   args.pop(idx);
   skipLines = args.pop(idx);
   #eprint("...");

keepMulti=False;
if "--keepMulti" in args:
   idx = args.index("--keepMulti");
   args.pop(idx);
   keepMulti = True;


writeTitleLine=False;
if "--writeTitleLine" in args:
   idx = args.index("--writeTitleLine");
   args.pop(idx);
   writeTitleLine = True;

if args[-1] == "-":
   args.pop(-1)
   inf = sys.stdin;
   isStream = True;
else:
   inf = args.pop(-1);
   isStream = False;
   inf = openSmart(inf,'r');


out = sys.stdout;

lnct = 0;
bodylnct = 0;
headerlnct=  0;

try:
   if writeTitleLine:
      out.write("CHROM\tSTART\tEND\tID\tREF\tALT\n");
   for line in inf:
      lnct = lnct + 1;
      if lnct % 1000 == 0:
         eprint(".", end='');
         if lnct % 5000 == 0:
            eprint(" ",end='');
         if lnct % 10000 == 0:
            cells = line.split("\t")
            eprint("["+str(lnct)+" lines complete] ["+time.strftime("%c")+"]")
      if not line.startswith("#"):
         bodylnct = bodylnct + 1;
         if bodylnct > skipLines:
            cells = line[:-1].split("\t")[0:5];
            pos = int(cells[1]);
            ref = cells[3];
            end = pos + len(ref) - 1;
            alt = cells[4]
            if not keepMulti:
               alt = cells[4].split(",")[0];
            out.write(cells[0]+"\t"+cells[1]+"\t"+str(end)+"\t"+cells[2]+"\t"+cells[3]+"\t"+alt+"\n");
      else:
         headerlnct = headerlnct + 1;

except IOError as e:
   if e.errno == errno.EPIPE:
      exit(0);

eprint("skipped "+str(headerlnct)+" header lines and "+str(skipLines)+" lines after that.");

out.close();

