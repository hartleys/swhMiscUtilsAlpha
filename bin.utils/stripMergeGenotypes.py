#!/usr/bin/env python
from __future__ import print_function

import sys,gzip,errno,time
#TODO-WRITE-HELPDOC

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

def indexOf(lst,x):
    try:
       y = lst.index(x);
       return y;
    except ValueError:
       return -1;

def isNumber(s):
    try:
        float(s)
        return True
    except ValueError:
        pass
    return False

eprint("STARTING");

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




eprint("STARTING");

args = sys.argv;

eprint("infile="+args[1]);
eprint("outfile="+args[2]);

convertToChr=False

if "--convertToChr" in args:
   args.pop(args.index("--convertToChr"));
   convertToChr=True
   


inf = openSmart(args[1],'r');
out = openSmart(args[2],'w');


lnct = 0
inHeader = True;
for line in inf:
   lnct = lnct + 1;
   if lnct % 1000 == 0:
      eprint(".", end='');
      if lnct % 5000 == 0:
         eprint(" ",end='');
      if lnct % 10000 == 0:
         eprint("["+str(lnct)+" lines complete] ["+time.strftime("%c")+"]")
   if line.startswith("##"):
      out.write(line);
   elif line.startswith("#"):
      out.write("##INFO=<ID=GT,Number=1,Type=STRING,Description=\"Genotype\",subType=\"GtStyle.unsplit\">\n");
      out.write("##INFO=<ID=DP,Number=1,Type=STRING,Description=\"Genotype depth\">\n");
      out.write("##INFO=<ID=GQ,Number=1,Type=STRING,Description=\"Genotype quality score\">\n");
      out.write("##INFO=<ID=AD,Number=R,Type=STRING,Description=\"Genotype quality score\">\n");
      #out.write("##INFO=<ID=FILTR,Number=1,Type=STRING,Description=\"GIAB quality filter field\">\n");
      cells = line[:-1].split("\t");
      out.write('\t'.join(cells[0:8])+"\n");
      inHeader = False;
      eprint("Finished with header...");
   else:
      cells = line[:-1].split("\t");
      fmtCells = cells[8].split(":");
      gtCells = cells[9].split(":");
      gt = gtCells[0];
      dpIdx = indexOf(fmtCells,"DP");
      gqIdx = indexOf(fmtCells,"GQ");
      adIdx = indexOf(fmtCells,"AD");
      if dpIdx != -1:
         dp = gtCells[dpIdx];
      else:
         dp = "0";
      if gqIdx != -1:
         gq = gtCells[gqIdx];
      else:
         gq = "0";
      if adIdx != -1:
         ad = gtCells[adIdx];
      else:
         ad = "0";
      linePrefix = ""
      if convertToChr:
         linePrefix="chr";
      
      out.write(linePrefix+'\t'.join(cells[0:8])+";GT="+gt+";DP="+dp+";GQ="+gq+";AD="+ad+"\n");

out.close();

eprint("DONE");



