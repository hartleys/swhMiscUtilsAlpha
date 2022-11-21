#!/usr/bin/env python
from __future__ import print_function

import sys,gzip,errno,time

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

eprint("STARTING");

args = sys.argv;

#+SHORTDESC="Converts VCF to tab-delimited table"
#+ HELPTEXT="vcf2table.py utility. Takes the info data from a VCF and makes a table. For testing purposes."$'\n'
#++HELPTEXT=" "$'\n'
#+ SYNTAXTEXT="vcf2table.py [options] myVcf.vcf.gz"$'\n'
#+ PARAMS="--help: displays syntax and help info."$'\n'
#++PARAMS="--man: synonym for --help."$'\n'
#++PARAMS="--infoColumnList <filename>: a comma-delimited list of INFO columns you want extracted, in order."$'\n'
#++PARAMS="--keepGT: Flag, keep genotype columns."$'\n'
#++PARAMS="--filterINFO [LIST]: ..."$'\n'
#++PARAMS="--filterGT [LIST]: ..."$'\n'
#++PARAMS="--tableGT: ..."$'\n'
#++PARAMS="--gtColumnList: ..."$'\n'
#++PARAMS="--emitOnlyNonref: ..."$'\n'
#++PARAMS="--infoIsNumericList [LIST]: Comma-delimited list of INFO tags that should be formatted as numerics, as opposed to quoted."$'\n'
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


def dequote(s):
    #If a string has single or double quotes around it, remove them.
    #Make sure the pair of quotes match.
    #If a matching pair of quotes is not found, return the string unchanged.
    if (s[0] == s[-1]) and s.startswith(("'", '"')):
        return s[1:-1]
    return s

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

infoCols = [];
presetInfoCols = False;
if "--keepGT" in args:
   keepGT = True;
else:
   keepGT = False;

filterInfoEQ = dict();
filterInfoNE = dict();
filterInfo = False;
if "--filterINFO" in args:
   idx = args.index("--filterINFO");
   filterListRawStrings = dequote(args.pop(idx+1)).split(",");
   args.pop(idx);
   filterInfo = True;
   filterPairsEQ  = [ x.split("==") for x in filterListRawStrings if len(x.split("==")) == 2]
   filterPairsNE  = [ x.split("!=") for x in filterListRawStrings if len(x.split("!=")) == 2]
   filterPairsERR = [ x for x in filterListRawStrings if len(x.split("==")) != 2 and len(x.split("!=")) != 2]
   if ( len(filterPairsERR) != 0 ):
      eprint("SYNTAX ERROR: filters must contain either == or !=");
   for fp in filterPairsEQ:
      filterInfoEQ[fp[0].lstrip()] = fp[1].lstrip();
      eprint("INFO Filtering for \""+fp[0].lstrip()+"\" EQ \""+fp[1].lstrip()+"\"");
   for fp in filterPairsNE:
      filterInfoNE[fp[0].lstrip()] = fp[1].lstrip();
      eprint("INFO Filtering for \""+fp[0].lstrip()+"\" NE \""+fp[1].lstrip()+"\"");

filterGenoEQ = dict();
filterGenoNE = dict();
filterGeno = False;
if "--filterGT" in args:
   idx = args.index("--filterGT");
   filterListRawStrings = dequote(args.pop(idx+1)).split(",");
   args.pop(idx);
   filterGeno = True
   keepGT = True
   filterPairsEQ  = [ x.split("==") for x in filterListRawStrings if len(x.split("==")) == 2]
   filterPairsNE  = [ x.split("!=") for x in filterListRawStrings if len(x.split("!=")) == 2]
   filterPairsERR = [ x.split("==") for x in filterListRawStrings if len(x.split("==")) != 2 and len(x.split("!=")) != 2]
   if ( len(filterPairsERR) != 0 ):
      eprint("SYNTAX ERROR: filters must contain either == or !=");
   for fp in filterPairsEQ:
      filterGenoEQ[fp[0].lstrip()] = fp[1].lstrip();
      eprint("Genotype Filtering for \""+fp[0].lstrip()+"\" EQ \""+fp[1].lstrip()+"\"");
   for fp in filterPairsNE:
      filterGenoNE[fp[0].lstrip()] = fp[1].lstrip();
      eprint("Genotype Filtering for \""+fp[0].lstrip()+"\" NE \""+fp[1].lstrip()+"\"");

tableGT=False;
if "--tableGT" in args:
   idx = args.index("--tableGT");
   args.pop(idx);
   tableGT = True;
   eprint("Making GT into a table...");

if "--infoColumnList" in args:
   idx = args.index("--infoColumnList");
   infoCols = args.pop(idx+1).split(",");
   args.pop(idx);
   presetInfoCols = True;

presetFmtCols = False;
fmtCols = [];

if "--gtColumnList" in args:
   idx = args.index("--gtColumnList");
   fmtCols = args.pop(idx+1).split(",");
   args.pop(idx);
   presetFmtCols = True;

emitOnlyNonref=False;
if "--emitOnlyNonref" in args:
   idx = args.index("--emitOnlyNonref");
   keepGT=True;
   emitOnlyNonref=True;
   filterGeno=True;
   args.pop(idx);

emitSampleIdsOnly=False;
if "--emitSampleIdsOnly" in args:
   idx = args.index("--emitSampleIdsOnly");
   emitSampleIdsOnly=True;
   args.pop(idx);

gtDelim="\t";
if "--gtDelim" in args:
   idx = args.index("--infoIsNumericList");
   gtDelim = args.pop(idx+1);
   args.pop(idx);

infoNumericSet = set();
if "--infoIsNumericList" in args:
   idx = args.index("--infoIsNumericList");
   args.pop(idx);
   infoNumericSet = set(args.pop(idx).split(","));


#vcfColNames = ["CHROM"  "POS"     "ID"      "REF"     "ALT"     "QUAL"    "FILTER"];
#vcfCols = [];
#presetInfoCols = False;
#if "--infoColumnList" in args:
#   idx = args.index("--infoColumnList");
#   vcfKeepColNames = args.pop(idx+1).split(",");
#   args.pop(idx);
#else:
#   vcfKeepColNames = vcfColNames;
#
#i = 0;
#for v in vcfColNames:
#   if v in vcfKeepColNames:
#      vcfCols = vcfCols + [i];
#   i = i + 1;

out = sys.stdout;

lnct = 0;
isOnHeader = True;
ALEN=16;

BLANKCHAR="."
NULLCHAR="-"
QUOTECHAR="\""
OPENQUOTE="=\""
LONGQUOTE="\""
CLOSEQUOTE="\"";
MISSINGNUMERICCHAR="-1";

infoFlagSet = set();

FORMATCOLUMN=8;

eprint("Starting Iteration ["+time.strftime("%c")+"]")

try:
   for line in inf:
      lnct = lnct + 1;
      if lnct % 1000 == 0 and not isOnHeader:
         eprint(".", end='');
         if lnct % 5000 == 0:
            eprint(" ",end='');
         if lnct % 10000 == 0:
            eprint("["+str(lnct)+" lines complete] ["+time.strftime("%c")+"]")
      if len(line) == 1:
         out.write("\n");
      elif line[0:2] != "##" and isOnHeader:
         #out.write("CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\t"+'\t'.join(infoCols));
         #if keepGT:
         #   out.write("\tFORMAT");
         #   if tableGT:
         #      out.write('\t'+'\t'.join(fmtCols));
         #out.write("\n");
         isOnHeader = False;
         eprint("Finished with header. ["+time.strftime("%c")+"]")
         titleLineCells = line[:-1].split("\t");
         sampleVarCounts = [0 for i in titleLineCells[(FORMATCOLUMN+1):]];
      elif line[0:7] == "##INFO=" :
         infoSplit = line[11:].split(",",2);
         infoNum = infoSplit[1][7:];
         infoKey = infoSplit[0];
         #eprint("header info line: key=\""+str(infoKey)+"\", num=\""+str(infoNum)+"\"");
         if not presetInfoCols:
            infoCols = infoCols + [infoKey];
         if infoNum == "0":
            infoFlagSet.add(infoKey);
            eprint("infoFlag: " + str(infoKey));
      elif line[0:9] == "##FORMAT=":
         infoSplit = line[13:].split(",",2);
         infoNum = infoSplit[1][7:];
         infoKey = infoSplit[0];
         #eprint("header info line: key=\""+str(infoKey)+"\", num=\""+str(infoNum)+"\"");
         if not presetFmtCols:
            fmtCols = fmtCols + [infoKey];
      elif line[0] != "#":
         cells = line[0:-1].split("\t");
         infocells = cells[7].split(";");
         infoPairs = [a.split('=',1) for a in infocells];
         infoNames = [a[0] for a in infoPairs];
         filterLine = False;
         if filterInfo:
            filterEqIdx = [ (a,infoNames.index(a)) for a in filterInfoEQ.keys() if a in infoNames]
            filterNeIdx = [ (a,infoNames.index(a)) for a in filterInfoNE.keys() if a in infoNames]
            for fp in filterEqIdx:
               if infoPairs[fp[1]][1] == filterInfoEQ[fp[0]]:
                 filterLine = True;
            for fp in filterNeIdx:
               if infoPairs[fp[1]][1] != filterInfoNE[fp[0]]:
                 filterLine = True;
         if not filterLine:
            #out.write('\t'.join(cells[0:7]));
            if not filterGeno:
               out.write( "\t"+ "\t".join(cells[9:]) );
            else:
               out.write("\t"+cells[FORMATCOLUMN]);
               fmt = cells[FORMATCOLUMN].split(":");
               filterEqIdx = [ (filterGenoEQ[a],fmt.index(a)) for a in filterGenoEQ.keys() if a in fmt]
               filterNeIdx = [ (filterGenoNE[a],fmt.index(a)) for a in filterGenoNE.keys() if a in fmt]
               #eprint("filterEqIdx: " + ",".join([str(a[0])+"/"+str(a[1]) for a in filterEqIdx]));
               #eprint("filterNeIdx: " + ",".join([str(a[0])+"/"+str(a[1]) for a in filterNeIdx]));
               if emitSampleIdsOnly:
                  out.write("\t");
               for idx in range(FORMATCOLUMN+1,len(cells)):
                  genoString = cells[idx];
                  filtGT = False;
                  geno = genoString.split(":");
                  for fp in filterEqIdx:
                     if len(geno) <= fp[1]:
                        if fp[0] == ".":
                           filtGT = True;
                     elif geno[fp[1]] == fp[0]:
                        filtGT = True;
                  for fp in filterNeIdx:
                     if len(geno) <= fp[1]:
                        if fp[0] != ".":
                           filtGT = True;
                     elif geno[fp[1]] != fp[0]:
                        filtGT = True;
                  if emitOnlyNonref:
                     if geno[0][0] == "1" or (len(geno[0]) > 2 and geno[0][2] == "1"):
                        filtGT = False;
                     else:
                        filtGT = True;
                  if not filtGT:
                     
except IOError as e:
   if e.errno == errno.EPIPE:
      exit(0);



out.close();

##INFO=<ID=

#
#  ../pullANN.py variants_annotated.subset.vcf.gz | less -S
#  ../pullANN.py variants_annotated.subset.vcf.gz | ../vcf2table.py - > test.txt
#  ../pullANN.py variants_annotated.subset.vcf.gz | ../vcf2table.py --infoColumnList ANN,SH_HGVSc_SHORT,SH_HGVSc_LONG,SH_HGVSc_SHORT_PC,SH_HGVSc_LONG_PC,SH_HGVSp_SHORT,SH_HGVSp_LONG - | gzip > test.txt.gz
#[C|sequence_feature|LOW|ITM2B|ENSG00000136156|topological_domain:Lumenal|ENST00000378565|protein_coding||c.565-192T>C||||||
#C|sequence_feature|LOW|ITM2B|ENSG00000136156|domain:BRICHOS|ENST00000378565|protein_coding||c.565-192T>C||||||
#C|sequence_feature|LOW|ITM2B|ENSG00000136156|disulfide_bond|ENST00000378565|protein_coding||c.565-192T>C||||||
#C|upstream_gene_variant|MODIFIER|ITM2B|ENSG00000136156|transcript|ENST00000463839|nonsense_mediated_decay||n.-2T>C|||||2532|
#C|intron_variant|MODIFIER|ITM2B|ENSG00000136156|transcript|ENST00000378565|protein_coding|4/5|c.565-192T>C||||||
#C|intron_variant|MODIFIER|ITM2B|ENSG00000136156|transcript|ENST00000378549|protein_coding|2/3|c.247-192T>C||||||
#C|intron_variant|MODIFIER|ITM2B|ENSG00000136156|transcript|ENST00000607866|nonsense_mediated_decay|4/4|n.*291-192T>C||||||]	NA

# ../vcf2table.py --infoColumnList SH_HGVSc_SHORT,SH_HGVSp_SHORT,1KG_AF annotated.v10.vcf.gz | less -S


