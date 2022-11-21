#!/usr/bin/env python
from __future__ import print_function

import sys,gzip,errno,time

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
#++PARAMS="--infoColumnList <tag1,tag2,...>: a comma-delimited list of INFO columns you want extracted, in order."$'\n'
#++PARAMS="--infoColumnReordering <tag1,tag2,...>: a comma-delimited list of INFO columns you want to be placed to the left of the other info columns, in order."$'\n'
#++PARAMS="--closeQuote: ..."$'\n'
#++PARAMS="--openQuote: ..."$'\n'
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

filterInfoGT = dict();
filterInfoLT = dict();

writeInfoFile=False
infoFile=""
if "--headerInfoFile" in args:
   idx = args.index("--headerInfoFile");
   infoFile=args.pop(idx+1);
   args.pop(idx);
   writeInfoFile=True;


noInfo = False;
fmt = False;
if "--noInfo" in args:
   idx = args.index("--noInfo");
   args.pop(idx);
   noInfo = True;

if "--fmt" in args:
   idx = args.index("--fmt");
   args.pop(idx);
   fmt = True;

if "--fmtOnly" in args:
   idx = args.index("--fmtOnly");
   args.pop(idx);
   fmt = True;
   noInfo = True;

initChar = ""
if "--commentOut" in args:
   idx = args.index("--commentOut");
   args.pop(idx);
   initChar = "#";
































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
   filterPairsGT  = [ x.split(">") for x in filterListRawStrings if len(x.split(">")) == 2]
   filterPairsLT  = [ x.split("<") for x in filterListRawStrings if len(x.split("<")) == 2]
   #filterPairsERR = [ x for x in filterListRawStrings if len(x.split("==")) != 2 and len(x.split("!=")) != 2]
   #if ( len(filterPairsERR) != 0 ):
   #   eprint("SYNTAX ERROR: filters must contain either == or !=");
   for fp in filterPairsEQ:
      filterInfoEQ[fp[0].strip()] = fp[1].strip();
      eprint("INFO Filtering out \""+fp[0].lstrip()+"\" EQ \""+fp[1].strip()+"\"");
   for fp in filterPairsNE:
      filterInfoNE[fp[0].strip()] = fp[1].strip();
      eprint("INFO Filtering out \""+fp[0].strip()+"\" NE \""+fp[1].strip()+"\"");
   for fp in filterPairsGT:
      filterInfoGT[fp[0].strip()] = float(fp[1].strip());
      eprint("INFO Filtering out \""+fp[0].strip()+"\" GT \""+fp[1].strip()+"\"");
   for fp in filterPairsLT:
      filterInfoLT[fp[0].strip()] = float(fp[1].lstrip());
      eprint("INFO Filtering out \""+fp[0].strip()+"\" LT \""+fp[1].strip()+"\"");

filterFilter = False;
if "--filterFilter" in args:
   idx = args.index("--filterFilter");
   args.pop(idx);
   filterFilter = True;

VERBOSE=False;
if "--verbose" in args:
   idx = args.index("--verbose");
   args.pop(idx);
   VERBOSE = True;
   eprint("--verbose reporting!");

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
   keepGT=True;
   eprint("Making GT into a table...");

gtLimit=-1;
if "--gtLimit" in args:
   idx = args.index("--gtLimit");
   gtLimit=int(args.pop(idx+1));
   args.pop(idx);
   eprint("Making GT into a table...");

if "--infoColumnList" in args:
   idx = args.index("--infoColumnList");
   infoCols = args.pop(idx+1).split(",");
   args.pop(idx);
   presetInfoCols = True;

infoOrdering = [];
if "--infoColumnReordering" in args:
   idx = args.index("--infoColumnReordering");
   infoOrdering = args.pop(idx+1).split(",");
   args.pop(idx);

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

gtTag="GT";
if "--gtTag" in args:
   idx = args.index("--gtTag");
   gtTag = args.pop(idx+1);
   args.pop(idx);

gtDelim="\t";
if "--gtDelim" in args:
   idx = args.index("--gtDelim");
   gtDelim = args.pop(idx+1);
   args.pop(idx);

infoNumericAuto = True;
infoNumericSet = set();
if "--infoIsNumericList" in args:
   idx = args.index("--infoIsNumericList");
   args.pop(idx);
   infoNumericSet = set(args.pop(idx).split(","));
   infoNumericAuto = False;

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
infolnct = 0;
fmtlnct = 0;
isOnHeader = True;
ALEN=16;

BLANKCHAR="."
NULLCHAR="-"
QUOTECHAR="\""
OPENQUOTE="=\""
#LONGQUOTE="\""
#CLOSEQUOTE="\"";
MISSINGNUMERICCHAR="-1";

if "--closeQuote" in args:
    idx = args.index("--closeQuote");
    args.pop(idx);
    QUOTECHAR = args.pop(idx)

if "--openQuote" in args:
    idx = args.index("--openQuote");
    args.pop(idx);
    OPENQUOTE = args.pop(idx)


infoFlagSet = set();

FORMATCOLUMN=8;
FILTERCOLUMN=6;

eprint("Starting Iteration ["+time.strftime("%c")+"]")

infoMetaData = dict();
fmtMetaData = dict();



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
      if( lnct < 10 and VERBOSE ):
         eprint("line["+lnct+"]:\""+line[:-1]+"\"")

      if len(line) == 1:
         out.write("\n");
         if( VERBOSE ):
            eprint("line["+lnct+"]: EMPTY")
      elif line[0:2] != "##" and isOnHeader:
         if( VERBOSE ):
            eprint("line["+lnct+"]: LAST HEADER LINE")
         isOnHeader = False;
         eprint("Finished with header. ["+time.strftime("%c")+"]")
         titleLineCells = line[:-1].split("\t");
         if len(infoOrdering) > 0:
            for infoTag in infoOrdering:
               if not infoTag in infoCols:
                  eprint("error: cannot find tag: "+infoTag);
               infoCols.remove(infoTag);
            infoCols = infoOrdering + infoCols;
         #out.write("CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\t"+'\t'.join(infoCols));
         #if keepGT:
         #   out.write("\tFORMAT");
         #   if tableGT:
         #      out.write('\t'+"sampleID"+'\t'+'\t'.join(fmtCols));
         #out.write("\n");
         if fmt and ( not noInfo) :
            out.write(initChar+"###--------------------------------------\n");
            out.write(initChar+"###INFO TAGS: \n")
         if  (not noInfo ) :
            out.write(initChar+"#infoTag\tnumber\ttype\tdescription\n");
            for infoTag in infoCols:
               md = infoMetaData[infoTag]
               out.write(initChar+infoTag+"\t"+md[1]+"\t"+md[2]+"\t"+md[3]+"\n");
         if fmt:
            if  not noInfo :
               out.write(initChar+"###--------------------------------------\n");
               out.write(initChar+"###FORMAT TAGS: \n")
            out.write(initChar+"#fmtTag\tnumber\ttype\tdescription\n");
            for fmtTag in fmtCols:
               md = fmtMetaData[fmtTag]
               out.write(initChar+fmtTag+"\t"+md[1]+"\t"+md[2]+"\t"+md[3]+"\n");
      elif line[0:7] == "##INFO=" :
         if( infolnct < 10 and VERBOSE ):
            eprint("line["+lnct+"]["+infolnct+"]: INFOLINE")
            infolnct = infolnct + 1;
         infoSplit = line[11:-1].split(",",3);
         infoTy  = infoSplit[2][5:];
         infoNum = infoSplit[1][7:];
         infoDesc = infoSplit[3][12:-1];
         infoKey = infoSplit[0];
         infoMetaData[infoKey] = [infoKey,infoNum,infoTy,infoDesc];
         #eprint("header info line: key=\""+str(infoKey)+"\", num=\""+str(infoNum)+"\"");
         if not presetInfoCols:
            infoCols = infoCols + [infoKey];
         if infoNum == "0":
            infoFlagSet.add(infoKey);
            eprint("infoFlag: " + str(infoKey));
         if infoNumericAuto and infoNum == "1" and (infoTy == "Integer" or infoTy == "Float"):
            infoNumericSet.add(infoKey);
      elif line[0:9] == "##FORMAT=":
         if( fmtlnct < 10 and VERBOSE ):
            eprint("line["+lnct+"]["+fmtlnct+"]: FMTLINE")
            infolnct = infolnct + 1;
         infoSplit = line[13:].split(",",3);
         infoTy  = infoSplit[2][5:];
         infoNum = infoSplit[1][7:];
         infoDesc = infoSplit[3][12:-1];
         infoKey = infoSplit[0];
         fmtMetaData[infoKey] = [infoKey,infoNum,infoTy,infoDesc];
         #eprint("header info line: key=\""+str(infoKey)+"\", num=\""+str(infoNum)+"\"");
         if not presetFmtCols:
            fmtCols = fmtCols + [infoKey];
      elif line[0] != "#":
         if( VERBOSE ):
             eprint("END OF HEADER! BREAK!")
         break;
except IOError as e:
   eprint("ERROR ERROR ERROR! open file failed!");
   eprint("ERROR ERROR ERROR: [\""+ "\",\"".join(e.args)+"\"]");
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
#C|intron_variant|MODIFIER|ITM2B|ENSG00000136156|transcript|ENST00000607866|nonsense_mediated_decay|4/4|n.*291-192T>C||||||]    NA

# ../vcf2table.py --infoColumnList SH_HGVSc_SHORT,SH_HGVSp_SHORT,1KG_AF annotated.v10.vcf.gz | less -S


