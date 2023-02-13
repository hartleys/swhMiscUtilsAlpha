#!/usr/bin/env python

from __future__ import print_function
#TODO: WRITE HELPDOC?

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
#args = ["acmg.table.2017-03-20.SWH.txt","Gene Symbol","ENSEMBL ID,Entrez Gene ID,Phenotype","./HC/BLAHBLAH","0","-"]
#./scripts/decodeColByName.py "acmg.table.2017-03-20.SWH.txt" "Gene Symbol" "ENSEMBL ID,Entrez Gene ID,Phenotype" "./scrap/test.in.txt" "SYMBOL" -
#./scripts/decodeColByName.py "acmg.table.2017-03-20.SWH.txt" "Gene Symbol" "ENSEMBL ID,Entrez Gene ID,Phenotype" "./HC/shAnnoAcmg/annotated.ACMG.wGRP.withGS.LLOF.table.txt" "SWH_SYMBOL" - | ./scripts/addQuotesToTable.py - - > ./HC/shAnnoAcmg/annotated.ACMG.wGRP.withGS.LLOF.info.table.txt



outfile = args.pop(-1);

fromCol2 = args.pop(-1);
if is_integer(fromCol2):
   fromColNum2 = int(fromCol2);
   colByNum2 = True;
else:
   fromColTitle2 = fromCol2;
   colByNum2 = False;

infile = args.pop(-1);

if "--listDelimiter" in args:
   idx = args.index("--listDelimiter");
   listDelim = args.pop(idx+1);
   args.pop(idx);
else:
   listDelim = ",";

if "--tableDelimiter" in args:
   idx = args.index("--tableDelimiter");
   tableDelimiter = args.pop(idx+1);
   args.pop(idx);
else:
   tableDelimiter = "\t";

if(is_integer(args[-2])):
   toColNums = [int(i) for i in args.pop(-1).split(listDelim)];
   fromColNum = int(args.pop(-1));
   toColTitles = ["X" + str(i) for i in range(0,len(toColNums))];
   colByNum = True;
else:
   colByNum = False;
   toColTitles = args.pop(-1).split(listDelim);
   fromColTitle = args.pop(-1);

decoderFile = args.pop(-1);

if "--singleHit" in args:
   assumeSingleHit = True;
else:
   assumeSingleHit = False;

#fromColTitle = ""
#
#
eprint("decoderFile  = \""+decoderFile+"\"");

if(colByNum):
   eprint("Reading decoder by column number...")
   eprint("fromColNum   = \""+str(fromColNum)+"\"");
   eprint("toColNums    = \""+str(toColNums)+"\"");
   eprint("toColTitles  = \""+str(toColTitles)+"\"");
else:
   eprint("Reading decoder by column title...")
   eprint("fromColTitle = \""+str(fromColTitle)+"\"");
   eprint("toColTitles  = \""+str(toColTitles)+"\"");

eprint("infile       = \""+infile+"\"");

if(colByNum2):
   eprint("Reading infile by column number...")
   eprint("fromColNum2 = \""+fromColNum2+"\"")
else:
   eprint("Reading infile by column title...")
   eprint("fromColTitle2 = \""+fromColTitle2+"\"")

eprint("outfile      = \""+outfile+"\"");

if decoderFile == "-":
   indataA = sys.stdin;
elif decoderFile[-3:] == ".gz":
   indataA = gzip.open(decoderFile,'r');
else:
   indataA = open(decoderFile,'r');

if infile == "-":
   indataB = sys.stdin;
elif infile[-3:] == ".gz":
   indataB = gzip.open(infile,'r');
else:
   indataB = open(infile,'r');

if outfile == "-":
   outdata = sys.stdout;
elif outfile[-3:] == ".gz":
   outdata = gzip.open(outfile,'w');
else:
   outdata = open(outfile,'w');


if not colByNum:
   titleLine = indataA.readline()[:-1].split(tableDelimiter);
   toColNums = [];
   for x in toColTitles:
      if not x in titleLine:
         raise Exception("ERROR: column title \"" + x + "\" not found in title line:\n"+str(titleLine));
      toColNums = toColNums + [titleLine.index(x)];
   if not fromColTitle in titleLine:
      raise Exception("ERROR: column title \"" + fromColTitle + "\" not found in title line:\n"+titleLine);
   fromColNum = titleLine.index(fromColTitle);
   eprint("fromColNum   = \""+str(fromColNum)+"\"");
   eprint("toColNums    = \""+str(toColNums)+"\"");

maxIdx = max(toColNums + [fromColNum])
lnct = 0;
fromDictA = collections.defaultdict(lambda: [set() for i in toColNums]);
countTableA = collections.defaultdict(lambda: 0);

for line in indataA:
   cells = line[:-1].split(tableDelimiter);
   lnct = lnct + 1;
   if len(cells) <= maxIdx:
      raise Exception("Bad decoder line: "+lnct+", less than "+str(maxIdx+1)+" columns");
   if cells[fromColNum] != "" and any([cells[i] != "" for i in toColNums]):
      for fr in cells[fromColNum].split(listDelim):
         countTableA[fr] = countTableA[fr] + 1;
         currList = fromDictA[fr];
         for i in range(0,len(toColNums)):
            currList[i] = currList[i] | set(cells[toColNums[i]].split(","));

eprint("Dict A:");
countTableTemp = collections.defaultdict(lambda: 0);
for k in sorted(countTableA.keys()):
   countTableTemp[countTableA[k]] = countTableTemp[countTableA[k]] + 1;

for ct in sorted(countTableTemp.keys()):
   eprint(str(ct)+"\t"+str(countTableTemp[ct]));

if not colByNum2:
   titleLine = indataB.readline()[:-1].split(tableDelimiter);
   if not fromColTitle2 in titleLine:
      raise Exception("ERROR: column title " + fromColTitle2 + " not found in title line:\n"+titleLine);
   fromColNum2 = titleLine.index(fromColTitle2);
   eprint("fromColNum2   = \""+str(fromColNum)+"\"");
   outdata.write(tableDelimiter.join(titleLine + toColTitles) + "\n");

for line in indataB:
   cells = line[:-1].split(tableDelimiter);
   outColStrings = ["" for i in toColNums];
   if cells[fromColNum2] != "":
      currSetList = [set() for i in toColNums];
      for fr in cells[fromColNum2].split(listDelim):
         dictOut = fromDictA[fr];
         for i in range(0,len(dictOut)):
            currSetList[i] = currSetList[i] | dictOut[i];
      outColStrings = [listDelim.join(sorted(currSet)) for currSet in currSetList];
   outdata.write(line[:-1] + tableDelimiter + tableDelimiter.join(outColStrings) + "\n");

outdata.close();



