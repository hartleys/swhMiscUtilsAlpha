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
#+ VERSION="0.0.5"


if len(args) == 0 or "--help" in args or "--man" in args:
   import inspect,commands,os
   SCRIPT_FILE_PATH=inspect.stack()[0][1];
   print(SCRIPT_FILE_PATH);
   SCRIPT_DIR_PATH=commands.getstatusoutput("dirname "+SCRIPT_FILE_PATH)[1]
   SCRIPT_FILENAME=commands.getstatusoutput("basename "+SCRIPT_FILE_PATH)[1]
   #print("SCRIPT_DIR_PATH="+SCRIPT_DIR_PATH);
   #print("SCRIPT_FILENAME="+SCRIPT_FILENAME);
   os.system(SCRIPT_DIR_PATH+"/helper/manForScript"+" "+SCRIPT_FILE_PATH)
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

out = sys.stdout;

lnct = 0;


FORMATCOLUMN=8;
FILTERCOLUMN=6;

infonames = [];

try:
   for line in inf:
      lnct = lnct + 1;
      if lnct % 1000 == 0 and not isOnHeader:
         eprint(".", end='');
         if lnct % 5000 == 0:
            eprint(" ",end='');
         if lnct % 10000 == 0:
            cells = line.split("\t")
            eprint("["+str(lnct)+" lines complete] ["+time.strftime("%c")+"] ["+cells[0]+":"+cells[1]+" "+cells[3]+">"+cells[4]+"]")
      cells = line[:-1].split("\t");
      if lnct == 1:
         out.write("##fileformat=VCFv4.1\n");
         infonames=cells;
         ##INFO=<ID=CLINVAR_ALLELEID,Number=.,Type=String,Description="Annotation tag ALLELEID from CLINVAR">
         for c in cells[8:]:
            out.write("##INFO=<ID="+c+",Number=.,Type=String,Description=\"Annotation tag ALLELEID from CLINVAR\">\n");
         out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");
      else:
         cleanCells = [ re.sub("[ ]+","_",re.sub("[;]",",",c)) for c in cells ];
         out.write(cleanCells[0:8]+"\t");
         for i in range(8,len(cells)):
            out.write(infonames[i]+"="+cleanCells[i]+";")
         out.write("\n");


out.close();


