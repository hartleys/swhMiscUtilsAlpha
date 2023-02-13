#!/usr/bin/env python
from __future__ import print_function

import sys,gzip,errno,time,re

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

def isNumber(s):
    try:
        float(s)
        return True
    except ValueError:
        pass
    return False

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
isOnHeader = True;
ALEN=16;
pathoidx=-1;

def cleanForVcf( s ):
   return re.sub('(^[_]+|[_]+$)','',re.sub('[^A-Za-z0-9_]+','',re.sub('[ :/\[\]|,;]+', '_', s)))

try:
   for line in inf:
      lnct = lnct + 1;
      if lnct % 1000 == 0:
         eprint(".", end='');
         if lnct % 5000 == 0:
            eprint(" ",end='');
         if lnct % 10000 == 0:
            cells = line.split("\t")
            eprint("["+str(lnct)+" lines complete] ["+time.strftime("%c")+"] ["+cells[0]+":"+cells[1]+" "+cells[3]+">"+cells[4]+"]")
      if lnct == 1:
         eprint("Writing header... ["+time.strftime("%c")+"]");
         #Chr	Start	End	Ref	Alt
         titleCells = line[:-1].split("\t");
         fieldTitleCells = titleCells[5:];
         out.write("##fileformat=VCFv4.1\n");
         for f in fieldTitleCells:
            out.write("##INFO=<ID="+cleanForVcf(f)+",Number=.,Type=String,Description=\"Field "+cleanForVcf(f)+" from Intervar\">\n");
         out.write("##INFO=<ID="+cleanForVcf("intervarCALL")+",Number=1,Type=String,Description=\"Intervar pathogenicity call\">\n");
         catLevels = ["PVS","PS","PM","PP","BA","BS","BP"];
         catLevelCts = [1,5,7,6,1,5,8]
         for ii in range(0,len(catLevels)):
            for jj in range(0,catLevelCts[ii]):
               out.write("##INFO=<ID="+cleanForVcf(catLevels[ii]+str(jj + 1))+",Number=1,Type=Integer,Description=\"ACMG Pathogenicity criterion "+catLevels[ii]+str(jj + 1)+" from Intervar\">\n");
         out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");
         eprint("done with header... ["+time.strftime("%c")+"]");
         #for ftc in fieldTitleCells:
         #   eprint("\""+ftc+"\"")
         pathoidx = fieldTitleCells.index(" InterVar: InterVar and Evidence ")
      else:
         cells = line[:-1].split("\t");
         fieldCells = cells[5:];
         out.write("chr"+cells[0]+"\t"+cells[1]+"\t.\t"+cells[3]+"\t"+cells[4]+"\t.\t.\t");
         for i in range(0,len(fieldTitleCells)):
            out.write(cleanForVcf(fieldTitleCells[i]) +"="+cleanForVcf(fieldCells[i])+";");


         pathString = cleanForVcf(fieldCells[pathoidx])
         pvsIdx = pathString.find("PVS1")
         out.write("intervarCALL" +"="+cleanForVcf(pathString[9:(pvsIdx-1)])+";");
         pathCells = pathString[pvsIdx:].split("_")
         for ii in range(0,len(catLevels)):
            lvl = catLevels[ii]
            if lvl != "PVS" and lvl != "BA":
               pathCells = pathCells[1:]
            elif lvl == "PVS":
               pathCells[0] = pathCells[0][4:]
            elif lvl == "BA":
               pathCells[0] = pathCells[0][3:]
            for jj in range(0,catLevelCts[ii]):
               out.write( catLevels[ii]+str(jj + 1)+"="+pathCells[0]+";");
               pathCells = pathCells[1:]


         out.write("\n");
except IOError as e:
   if e.errno == errno.EPIPE:
      exit(0);

out.close();

