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
#++PARAMS="--emitOnlyNonref: for use with tableGT, with this option only samples with at least one '1' allele will have a line. This is sometimes referred to as a 'longTable'"$'\n'
#++PARAMS="--emitAnyNonref: For use with tableGT, with this option only samples with any non-reference allele will have a line."$'\n'
#++PARAMS="--emitSampleIdsOnly: ..."$'\n'
#++PARAMS="--emitOnlySamplesOnList [LIST]: For use with tableGT, with this option only samples on the supplied list will have a line."$'\n'
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

printInfoFile=False
if "--header" in args or "-H" in args:
   if "--header" in args:
      idx = args.index("--header");
   else:
      idx = args.index("-H");
   args.pop(idx);
   printInfoFile=True;

urlCols = dict();
while "--urlColumn" in args:
   idx = args.index("--urlColumn");
   currUrlCol=args.pop(idx+1);
   args.pop(idx);
   currUrlCells = currUrlCol.split("|")
   urlCols[currUrlCells[0]] = currUrlCells[1:];
   eprint("adding URL: \""+currUrlCells[0]+"\"[ \""+"\",\"".join(currUrlCells[1:])+"\"]");

eprint("urlCols:"+str(urlCols));

writeExtendedInfoFile=False
if "--extendedHeaderInfoFile" in args:
   idx = args.index("--extendedHeaderInfoFile");
   infoFile=args.pop(idx+1);
   args.pop(idx);
   writeExtendedInfoFile=True;

extendedInfoFileExampleCt=10;

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

flatTableGT=False;
if "--flatTableGT" in args:
   idx = args.index("--flatTableGT");
   args.pop(idx);
   flatTableGT = True;
   keepGT=True;
   eprint("Making GT into a FLAT table...");

gtLimit=-1;
if "--gtLimit" in args:
   idx = args.index("--gtLimit");
   gtLimit=int(args.pop(idx+1));
   args.pop(idx);
   eprint("Making GT into a table...");


truncateFields=-1;
if "--truncateFields" in args:
   idx = args.index("--truncateFields");
   truncateFields=int(args.pop(idx+1));
   args.pop(idx);
   eprint("Truncating fields to length: "+str(truncateFields));


if "--infoColumnList" in args:
   idx = args.index("--infoColumnList");
   infoCols = args.pop(idx+1).split(",");
   args.pop(idx);
   presetInfoCols = True;


keepsamplist = [];
if "--keepSamples" in args:
   idx = args.index("--keepSamples");
   args.pop(idx)
   keepsamplist = args.pop(idx).split(",");


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

flatBufferLength=-1;
if "--flatBufferLength" in args:
   idx = args.index("--flatBufferLength");
   flatBufferLength = int(args.pop(idx+1));
   args.pop(idx);

emitOnlyNonref=False;
if "--emitOnlyNonref" in args:
   idx = args.index("--emitOnlyNonref");
   keepGT=True;
   emitOnlyNonref=True;
   filterGeno=True;
   args.pop(idx);

emitPctRef=-1;
if "--emitPctRef" in args:
   idx = args.index("--emitPctRef");
   args.pop(idx);
   emitPctRef=float(args.pop(idx));
   eprint("emitPctRef="+str(emitPctRef));

emitOnlySamplesOnList=False;
if "--emitOnlySamplesOnList" in args:
   idx = args.index("--emitOnlySamplesOnList");
   keepGT=True;
   filterGeno=True;
   args.pop(idx);
   emitOnlySamplesOnList=args.pop(idx).split(",");

eprint("TEST emitPctRef="+str(emitPctRef));


emitAnyNonref=False;
if "--emitAnyNonref" in args:
   idx = args.index("--emitAnyNonref");
   keepGT=True;
   emitAnyNonref=True;
   filterGeno=True;
   args.pop(idx);

emitOnlySimpleHet=False;
if "--emitOnlySimpleHet" in args:
   idx = args.index("--emitOnlySimpleHet");
   keepGT=True;
   emitOnlySimpleHet=True;
   filterGeno=True;
   args.pop(idx);


emitOnlyAny=False;
if "--emitOnlyAny" in args:
   idx = args.index("--emitOnlyAny");
   args.pop(idx);
   emitOnlyAny= args.pop(idx).split(",")
   for ee in range(0,len(emitOnlyAny)):
      emitOnlyAny[ee] = emitOnlyAny[ee].split(":");
      if len( emitOnlyAny[ee] ) != 2:
          eprint("ERROR: each comma-delimited element in emitOnly must be of form tag:style where tag is a Genotype field and style is either AnyAlt, SimpleHet, AnyNonref, or 1");
          exit(1)
      if emitOnlyAny[ee][1] not in ["AnyAlt","SimpleHet","AnyNonref","1","0","NotAnyAlt"]:
          if emitOnlyAny[ee][1].startswith("eq."):
             emitOnlyAny[ee][1] = emitOnlyAny[ee][1][3:];
          else:
             eprint("ERROR: each comma-delimited element in emitOnly must be of form tag:style where tag is a Genotype field and style is either AnyAlt, SimpleHet, AnyNonref, or 1");
             exit(1)



emitOnlyAll=False;
if "--emitOnlyAll" in args:
   idx = args.index("--emitOnlyAll");
   args.pop(idx);
   emitOnlyAll= args.pop(idx).split(",")
   for ee in range(0,len(emitOnlyAll)):
      emitOnlyAll[ee] = emitOnlyAll[ee].split(":");
      if len( emitOnlyAll[ee] ) != 2:
          eprint("ERROR: each comma-delimited element in emitOnly must be of form tag:style where tag is a Genotype field and style is either AnyAlt, SimpleHet, AnyNonref, or 1");
          exit(1)
      if emitOnlyAll[ee][1] not in ["AnyAlt","SimpleHet","AnyNonref","1","0","NotAnyAlt"]:
          if emitOnlyAll[ee][1].startswith("eq."):
             emitOnlyAll[ee][1] = emitOnlyAll[ee][1][3:];
          else:
             eprint("ERROR: each comma-delimited element in emitOnly must be of form tag:style where tag is a Genotype field and style is either AnyAlt, SimpleHet, AnyNonref, or 1");
             exit(1)

if emitOnlyAny:
   eprint("EMIT ONLY ANY...")
   for (gtt,style) in emitOnlyAny:
       eprint("EMIT ONLY ANY: "+gtt+"  :  "+style);

if emitOnlyAll:
   eprint("EMIT ONLY ALL...")
   for (gtt,style) in emitOnlyAll:
       eprint("EMIT ONLY ALL: "+gtt+"  :  "+style);



emitSampleIdsOnly=False;
if "--emitSampleIdsOnly" in args:
   idx = args.index("--emitSampleIdsOnly");
   emitSampleIdsOnly=True;
   args.pop(idx);

gtTag="GT";
if "--gtTag" in args:
   idx = args.index("--gtTag");
   gtTag = args.pop(idx+1).split(",");
   args.pop(idx);

replaceUnderscores=False;
if "--replaceUnderscores" in args:
   idx = args.index("--replaceUnderscores");
   args.pop(idx);
   replaceUnderscores = True;

replaceInfoUnderscores=False;
if "--replaceInfoUnderscores" in args:
   idx = args.index("--replaceInfoUnderscores");
   args.pop(idx);
   replaceInfoUnderscores = True;


DEFAULT_NUMASNUMS = True;

numsAsNums=DEFAULT_NUMASNUMS;

#numsAsNums=False;
if "--numsAsNums" in args:
   idx = args.index("--numsAsNums");
   args.pop(idx);
   numsAsNums = True;

quoteNums=False;
if "--quoteNumbers" in args:
   idx = args.index("--quoteNumbers");
   args.pop(idx);
   numsAsNums = False;
   quoteNums = True;

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

if "--infoIsNumericSet" in args:
   idx = args.index("--infoIsNumericSet");
   args.pop(idx);
   infoNumericSet = set(args.pop(idx).split(","));


annFields = [];
while "--snpEffAnnField" in args:
   idx = args.index("--snpEffAnnField");
   args.pop(idx);
   annFields = annFields + [args.pop(idx)];

effFields = [];
while "--snpEffEffField" in args:
   idx = args.index("--snpEffEffField");
   args.pop(idx);
   effFields = effFields + [args.pop(idx)];


for annID in annFields:
   eprint("annID="+annID);

hasSampleDecoder = False
sampDict = collections.defaultdict( list )
sampDictTitleCells = [];
if "--sampleDecoder" in args:
   idx = args.index("--sampleDecoder");
   sampDictFile=args.pop(idx+1);
   args.pop(idx);
   lnctxx = 0;
   for line in openSmart(sampDictFile,'r'):
      cells = line[:-1].split("\t");
      if lnctxx == 0:
          sampDictTitleCells = cells[1:];
      else:
          sampDict[cells[0]] = cells[1:];
      lnctxx += 1
   hasSampleDecoder = True;

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
#LONGQUOTE="\""
#CLOSEQUOTE="\"";
MISSINGNUMERICCHAR="-1";

if "--missingVal" in args:
   idx = args.index("--missingVal");
   args.pop(idx);
   NULLCHAR=args.pop(idx)
   MISSINGNUMERICCHAR=NULLCHAR


if "--closeQuote" in args:
    idx = args.index("--closeQuote");
    args.pop(idx);
    QUOTECHAR = args.pop(idx)

if "--openQuote" in args:
    idx = args.index("--openQuote");
    args.pop(idx);
    OPENQUOTE = args.pop(idx)

if "--noQuote" in args:
    idx = args.index("--noQuote");
    args.pop(idx);
    OPENQUOTE = ""
    QUOTECHAR= ""


infoFlagSet = set();
infoFlagWarnSet = set();

FORMATCOLUMN=8;
FILTERCOLUMN=6;

eprint("Starting Iteration ["+time.strftime("%c")+"]")

infoMetaData = dict();

annFieldSubCols = ["ALT","EFF","LVL","geneSymbol","geneID","featureType","featureID","txType","rankTotal","HGVSc","HGVSp","cPos","cdsPos","protPos","dist","warn"]
effFieldSubCols = ["EFF","Impact","LVL","funcClass","codonChange","aaChange","aaLen","geneName","txType","geneCoding","txID","exonRank","genotype","err","warn"]

##Effect ( Effect_Impact | Functional_Class | Codon_Change | Amino_Acid_Change| Amino_Acid_length | Gene_Name | Transcript_BioType | Gene_Coding | Transcript_ID | Exon_Rank  | Genotype [ | ERRORS | WARNINGS ] )

#if replaceInfoUnderscores:
#  def prepInfo(s):
#     return re.sub("[_]"," ",s);
#else:
#  def prepInfo(s):
#     return s
#

WARN_GENO_MISSING = 0;
WARN_GENO_MISSING_LIMIT = 10;

try:
   for line in inf:
      if len(line) == 1:
         out.write("\n");
      elif line[0:2] != "##" and isOnHeader:
         isOnHeader = False;
         eprint("Finished with header. ["+time.strftime("%c")+"]")
         titleLineCells = line[:-1].split("\t");
         if len(infoOrdering) > 0:
            for infoTag in infoOrdering:
               if not infoTag in infoCols:
                  eprint("error: cannot find tag: "+infoTag);
               infoCols.remove(infoTag);
            infoCols = infoOrdering + infoCols;
         infoColStrings = infoCols;
         if replaceUnderscores:
            infoColStrings = [x.replace("_"," ") for x in infoCols ]
         out.write("CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\t"+'\t'.join(infoColStrings));
         for annID in annFields:
            out.write("\t"+"\t".join([ annID+"_"+ annSub for annSub in annFieldSubCols ] ));
         for annID in effFields:
            out.write("\t"+"\t".join([ annID+"_"+ annSub for annSub in effFieldSubCols ] ));
         if keepGT:
            if tableGT :
               out.write('\tIDX\t'+"sampleID"+'\t'+'\t'.join(fmtCols + sampDictTitleCells));
            elif flatTableGT:
               if emitOnlySimpleHet or emitOnlyNonref or emitAnyNonref or filterGeno or emitOnlySamplesOnList or emitOnly:
                  out.write("\tFORMAT");
                  if flatBufferLength == -1:
                     for idx in range(FORMATCOLUMN+1,len(titleLineCells)):
                        #sampid = str(idx - FORMATCOLUMN)
                        #if len(keepsamplist) == 0 or sampid in keepsamplist:
                        out.write('\t'+'\t'.join( ["S_SAMPID"]+ ["S_" + str(idx - FORMATCOLUMN) + "_" + f for f in fmtCols] ));
                  else:
                     for idx in range(1,flatBufferLength + 1):
                        #sampid = str(idx - FORMATCOLUMN)
                        #if len(keepsamplist) == 0 or sampid in keepsamplist:
                        out.write('\t'+'\t'.join( ["S_SAMPID"]+["S_" + str(idx) + "_" + f for f in fmtCols] ));
               else :
                  out.write("\tFORMAT");
                  for idx in range(FORMATCOLUMN+1,len(titleLineCells)):
                     sampid = titleLineCells[idx]
                     if len(keepsamplist) == 0 or sampid in keepsamplist:
                        out.write('\t'+'\t'.join([titleLineCells[idx] + "_" + f for f in fmtCols]));
         out.write("\n");

         if writeInfoFile:
            infoWriter = open(infoFile,'w');
            infoWriter.write("infoTag\tnumber\ttype\tdescription\n");
            for infoTag in infoCols:
               md = infoMetaData[infoTag]
               infoWriter.write(infoTag+"\t"+md[1]+"\t"+md[2]+"\t"+md[3]+"\n");
            infoWriter.close();
         elif printInfoFile:
            out.write("infoTag\tnumber\ttype\tdescription\n");
            for infoTag in infoCols:
               md = infoMetaData[infoTag]
               out.write(infoTag+"\t"+md[1]+"\t"+md[2]+"\t"+md[3]+"\n");
            out.close();
            sys.exit();
      elif line[0:7] == "##INFO=" :
         infoSplit = line[11:-1].split(",",3);
         infoTy  = infoSplit[2][5:];
         infoNum = infoSplit[1][7:];
         infoDesc = infoSplit[3][12:-1];
         infoKey = infoSplit[0];
         #infoMetaData[infoKey] = [infoKey,infoNum,infoTy,infoDesc];
         if infoTy == "Float":
           md = [infoKey,infoNum,infoTy,infoDesc,0, float('inf'),float('-inf'), float('inf'),float('-inf'), set()];
         elif infoTy == "String":
           md = [infoKey,infoNum,infoTy,infoDesc,0, float('inf'),float('-inf'), 0,0, set()];
         elif infoTy == "Integer":
           md = [infoKey,infoNum,infoTy,infoDesc,0, float('inf'),float('-inf'), float('inf'),float('-inf'), set()];
         else:
           md = [infoKey,infoNum,infoTy,infoDesc,0, float('inf'),float('-inf'), float('inf'),float('-inf'), set()];
         infoMetaData[infoKey] = md;
         #md = [infoKey,infoNum,infoTy,infoDesc,missingCt,lenMin,lenMax,min,max,examples]
         #eprint("header info line: key=\""+str(infoKey)+"\", num=\""+str(infoNum)+"\"");
         if not presetInfoCols:
            infoCols = infoCols + [infoKey];
         if infoNum == "0" :
            #and ( not (infoKey in infoNumericSet)) :
            infoFlagSet.add(infoKey);
            eprint("infoFlag: " + str(infoKey));
         if infoNumericAuto and infoNum == "1" and (infoTy == "Integer" or infoTy == "Float"):
            infoNumericSet.add(infoKey);
      elif line[0:9] == "##FORMAT=":
         infoSplit = line[13:].split(",",2);
         infoNum = infoSplit[1][7:];
         infoKey = infoSplit[0];
         #eprint("header info line: key=\""+str(infoKey)+"\", num=\""+str(infoNum)+"\"");
         if not presetFmtCols:
            fmtCols = fmtCols + [infoKey];
      elif line[0] != "#":
         lnct = lnct + 1;
         if lnct < 10:
           eprint(".",end='');
         elif lnct == 10:
           eprint(" done 10 lines");
         elif lnct < 100 and lnct % 10 == 0:
           eprint(".",end='');
         elif lnct == 100:
           eprint(" done 100 lines");
           eprint(" ", end='');
         elif lnct % 1000 == 0:
           eprint(".", end='');
           if lnct % 5000 == 0:
              eprint(" ",end='');
           if lnct % 10000 == 0:
              cells = line.split("\t")
              eprint("["+str(lnct)+" lines complete] ["+time.strftime("%c")+"] ["+cells[0]+":"+cells[1]+" "+cells[3]+">"+cells[4]+"]")
         cells = line[0:-1].split("\t");
         infocells = cells[7].split(";");
         infoPairs = [a.split('=',1) for a in infocells];
         infoNames = [a[0] for a in infoPairs];
         filterLine = False;
         if filterFilter and ( cells[FILTERCOLUMN] != "PASS" and cells[FILTERCOLUMN] != "." ):
            filterLine = True;
         if filterInfo:
            filterEqIdx = [ (a,infoNames.index(a)) for a in filterInfoEQ.keys() if a in infoNames]
            filterNeIdx = [ (a,infoNames.index(a)) for a in filterInfoNE.keys() if a in infoNames]
            filterGtIdx = [ (a,infoNames.index(a)) for a in filterInfoGT.keys() if a in infoNames]
            filterLtIdx = [ (a,infoNames.index(a)) for a in filterInfoLT.keys() if a in infoNames]

            for fp in filterEqIdx:
               if infoPairs[fp[1]][1] == filterInfoEQ[fp[0]]:
                 filterLine = True;
            for fp in filterNeIdx:
               if infoPairs[fp[1]][1] != filterInfoNE[fp[0]]:
                 filterLine = True;
            for fp in filterGtIdx:
               if isNumber(infoPairs[fp[1]][1]):
                  if float(infoPairs[fp[1]][1]) > filterInfoGT[fp[0]]:
                     filterLine = True;
            for fp in filterLtIdx:
               if isNumber(infoPairs[fp[1]][1]):
                  if float(infoPairs[fp[1]][1]) < filterInfoLT[fp[0]]:
                     filterLine = True;
         if not filterLine:
            outline = '\t'.join(cells[0:7]);
            #out.write('\t'.join(cells[0:7]));
            #infoMetaData[infoKey] = [0 infoKey,1 infoNum,2 infoTy,3 infoDesc,4 missingCt,5 lenMin,6 lenMax,7 min,8 max,9 examples]
            for c in infoCols:

               if c.startswith("IDX:"):
                  cc = c.split(":")[2];
                  if cc in infoNames:
                     iix = intOrZero(c.split(":")[1]);
                     idx = infoNames.index(cc)
                     tagValue = infoPairs[idx][1]
                     tagCells = tagValue.split(",");
                     if iix < len(tagCells):
                        outline += "\t"+OPENQUOTE+tagCells[iix]+QUOTECHAR
                     else :
                        outline += "\t"+OPENQUOTE+NULLCHAR+QUOTECHAR;
                  else:
                     outline += "\t"+OPENQUOTE+NULLCHAR+QUOTECHAR;
               elif c.startswith("URL:"):
                  cc = c.split(":")[1];
                  if cc in infoNames:
                     idx = infoNames.index(cc)
                     urlPrefix = urlCols[cc][0];
                     urlSuffix = urlCols[cc][1];
                     tagValue = infoPairs[idx][1]
                     if len(urlCols[cc]) == 3:
                        urlInfix = urlCols[cc][2];
                        tagCells = tagValue.split(",");
                        outline += "\t=HYPERLINK(\""+urlPrefix+urlInfix.join(tagCells)+urlSuffix+"\",\""+tagValue+"\")"
                     else:
                        outline += "\t=HYPERLINK(\""+urlPrefix+tagValue+urlSuffix+"\",\""+tagValue+"\")"
                  else:
                     outline += "\t"+OPENQUOTE+NULLCHAR+QUOTECHAR;
               elif c in infoNames:
                  if c in infoFlagSet:
                     idx = infoNames.index(c)
                     if len(infoPairs[idx]) == 1:
                        outline += "\t"+OPENQUOTE+"TRUE"+QUOTECHAR
                        #out.write("\t"+OPENQUOTE+"TRUE"+QUOTECHAR);
                     else:
                        if ( c not in infoFlagWarnSet ):
                           infoFlagWarnSet.add(c)
                           eprint("\n   WARNING: VCF FORMAT SPEC VIOLATION: INFO tag \""+c+"\" is listed as a flag but actually appears to have an assigned value for at least one line!");
                           eprint("        Continuing with processing, but this column may have a mix of logic and nonlogic values.")
                           eprint("        Depending on what this INFO field is supposed to mean, this may or may not be the desired behavior!")
                           eprint("        tag: \""+c+"\" = \""+infoPairs[idx][1]+"\"")
                        if (numsAsNums and isNumber(infoPairs[idx][1])):
                           outline += "\t"+infoPairs[idx][1]
                        else:
                           outline += "\t"+OPENQUOTE+infoPairs[idx][1]+QUOTECHAR
                        #out.write("\t"+OPENQUOTE+infoPairs[idx][1]+QUOTECHAR);
                  elif c in infoNumericSet:
                     idx = infoNames.index(c)
                     if len(infoPairs[idx]) == 1:
                        outline += "\t"+MISSINGNUMERICCHAR
                        #out.write("\t"+MISSINGNUMERICCHAR);
                     elif infoPairs[idx][1] == ".":
                        outline += "\t"+MISSINGNUMERICCHAR
                        #out.write("\t"+MISSINGNUMERICCHAR);
                     else:
                        outline += "\t"+infoPairs[idx][1]
                        #out.write("\t"+infoPairs[idx][1]);
                  else:
                     idx = infoNames.index(c)
                     if len(infoPairs[idx]) == 1:
                        outline += "\t"+OPENQUOTE+BLANKCHAR+QUOTECHAR
                        #out.write("\t"+OPENQUOTE+BLANKCHAR+QUOTECHAR);
                     else:
                        if (numsAsNums and isNumber(infoPairs[idx][1])):
                           outline += "\t"+infoPairs[idx][1]
                        elif truncateFields == -1:
                           outline += "\t"+OPENQUOTE+infoPairs[idx][1]+QUOTECHAR
                        else:
                           outline += "\t"+OPENQUOTE+infoPairs[idx][1][0:truncateFields]+QUOTECHAR
                        #outline += "\t"+OPENQUOTE+prepInfo(infoPairs[idx][1])+QUOTECHAR
                        #out.write("\t"+OPENQUOTE+infoPairs[idx][1]+QUOTECHAR);
                  if writeExtendedInfoFile:
                     infoNum = infoMetaData[c][1]
                     infoTy = infoMetaData[c][2]
                     if len(infoPairs[idx]) > 1:
                        infoList = infoPairs[idx][1].split(",");
                        if infoNum != "1":
                           listLen = len(infoList);
                           infoMetaData[c][5] = min(infoMetaData[c][5],listLen);
                           infoMetaData[c][6] = max(infoMetaData[c][6],listLen);
                        if infoTy == "Float":
                           floatValues = [floatOrZero(x) for x in infoList]
                           infoMetaData[c][7] = min(infoMetaData[c][7],min(floatValues));
                           infoMetaData[c][8] = max(infoMetaData[c][8],max(floatValues));
                        elif infoTy == "Integer":
                           intValues = [intOrZero(x) for x in infoList]
                           infoMetaData[c][7] = min(infoMetaData[c][7],min(intValues));
                           infoMetaData[c][8] = max(infoMetaData[c][8],max(intValues));
                        else:
                           if (infoList[0] not in infoMetaData[c][9]):
                              if len(infoMetaData[c][9]) <= extendedInfoFileExampleCt :
                                 infoMetaData[c][9].add(infoList[0])
                              elif random.randint(0,lnct) == 1:
                                 infoMetaData[c][9].pop();
                                 infoMetaData[c][9].add(infoList[0]);
               else:
                  infoMetaData[c][4] += 1;
                  if c in infoFlagSet:
                     outline += "\t"+OPENQUOTE+"FALSE"+QUOTECHAR
                     #out.write("\t"+OPENQUOTE+"FALSE"+QUOTECHAR);
                  elif c in infoNumericSet:
                     outline += "\t"+MISSINGNUMERICCHAR
                     #out.write("\t"+MISSINGNUMERICCHAR);
                  else:
                     outline += "\t"+OPENQUOTE+NULLCHAR+QUOTECHAR
                     #out.write("\t"+OPENQUOTE+NULLCHAR+QUOTECHAR);
            for annID in annFields:
               if annID in infoNames:
                  idx = infoNames.index(annID)
                  #re.split(pattern, string)
                  annData = infoPairs[idx][1].split(",")[0].split("|")
                  outline += "\t"+"\t".join( [ OPENQUOTE+xx+QUOTECHAR for xx in annData ] );
               else:
                  for zz in range(0,len(annFieldSubCols)):
                     outline += "\t"+OPENQUOTE+NULLCHAR+QUOTECHAR
            for annID in effFields:
               if annID in infoNames:
                  idx = infoNames.index(annID)
                  annData = re.split("[|()]", infoPairs[idx][1].split(",")[0])
                  annIdx = 0;
                  for xx in annData:
                     if annIdx < len(effFieldSubCols):
                        annIdx = annIdx + 1;
                        if xx == "":
                           outline += "\t"+OPENQUOTE+"."+QUOTECHAR;
                        else:
                           outline += "\t"+OPENQUOTE+xx+QUOTECHAR;
                  while annIdx < len(effFieldSubCols):
                     annIdx = annIdx + 1;
                     outline += "\t"+OPENQUOTE+NULLCHAR+QUOTECHAR;
               else:
                  for zz in range(0,len(annFieldSubCols)):
                     outline += "\t"+OPENQUOTE+NULLCHAR+QUOTECHAR
            if(keepGT):

               if (not filterGeno) and (not emitSampleIdsOnly) and (not tableGT) and (not flatTableGT):
                  out.write( outline+"\t"+cells[FORMATCOLUMN]+"\t"+ "\t".join(cells[9:]) );
               else:
                  if not tableGT:
                     out.write(outline+"\t"+cells[FORMATCOLUMN]);
                  fmt = cells[FORMATCOLUMN].split(":");

                  if filterGeno:
                     filterEqIdx = [ (filterGenoEQ[a],fmt.index(a)) for a in filterGenoEQ.keys() if a in fmt]
                     filterNeIdx = [ (filterGenoNE[a],fmt.index(a)) for a in filterGenoNE.keys() if a in fmt]
                  else:
                     filterEqIdx = [];
                     filterNeIdx = [];

                  #if gtTag in fmt:
                  if True:
                     #gtTagIdx = fmt.index(gtTag);
                     printGtCt = 0;
                     if emitSampleIdsOnly:
                        out.write("\t");
                     for idx in range(FORMATCOLUMN+1,len(cells)):
                        genoString = cells[idx];
                        filtGT = False;
                        geno = genoString.split(":");
                        sampid = titleLineCells[idx]
                        if len(keepsamplist) > 0:
                           if sampid not in keepsamplist:
                              filtGT = True;
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
                              #emitOnlySimpleHet
                        if emitOnlyAny:
                           isAlt=False;
                           for (gtt,gtStyle) in emitOnlyAny:
                              if gtt in fmt:
                                 gtTagIdx = fmt.index(gtt);
                                 if gtStyle == "SimpleHet" and ( geno[gtTagIdx] == "0/1" or geno[gtTagIdx] == "1/0" ):
                                    isAlt = True;
                                 elif gtStyle == "AnyAlt" and ( "1" in re.split("[|/]",geno[gtTagIdx]) ):
                                    isAlt = True;
                                 elif gtStyle == "NotAnyAlt" and ( "1" not in re.split("[|/]",geno[gtTagIdx]) ):
                                    currAlt = True;
                                 elif gtStyle == "AnyNonref":
                                    genoCells = re.split("[|/]",geno[gtTagIdx])
                                    for gg in genoCells:
                                       if gg not in [".","0"]:
                                          isAlt = True;
                                 elif gtStyle == "1" and ( geno[gtTagIdx] == "1" ):
                                    isAlt = True;
                                 elif gtStyle == "0" and ( geno[gtTagIdx] == "0" ):
                                    currAlt = True;
                                 elif geno[gtTagIdx] == gtStyle:
                                    currAlt = True;
                           anyObsAlt=isAlt
                           if ( not anyObsAlt ) and emitPctRef < 0:
                              filtGT = True;
                           elif (not anyObsAlt) and random.random() > emitPctRef:
                              filtGT = True;
                        if emitOnlyAll:
                           isAlt=True;
                           for (gtt,gtStyle) in emitOnlyAll:
                              if gtt in fmt:
                                 gtTagIdx = fmt.index(gtt);
                                 currAlt=False;
                                 if gtStyle == "AnyAlt" and ( geno[gtTagIdx] == "0/1" or geno[gtTagIdx] == "1/0" ):
                                    currAlt = True;
                                 elif gtStyle == "AnyAlt" and ( "1" in re.split("[|/]",geno[gtTagIdx]) ):
                                    currAlt = True;
                                 elif gtStyle == "NotAnyAlt" and ( "1" not in re.split("[|/]",geno[gtTagIdx]) ):
                                    currAlt = True;
                                 elif gtStyle == "AnyNonref":
                                    genoCells = re.split("[|/]",geno[gtTagIdx])
                                    for gg in genoCells:
                                       if gg not in [".","0"]:
                                          currAlt = True;
                                 elif gtStyle == "1" and ( geno[gtTagIdx] == "1" ):
                                    currAlt = True;
                                 elif gtStyle == "0" and ( geno[gtTagIdx] == "0" ):
                                    currAlt = True;
                                 elif geno[gtTagIdx] == gtStyle:
                                    currAlt = True;
                                 isAlt = isAlt and currAlt;
                              else:
                                 isAlt=False;
                           anyObsAlt=isAlt
                           if ( not anyObsAlt ) and emitPctRef < 0:
                              filtGT = True;
                           elif (not anyObsAlt) and random.random() > emitPctRef:
                              filtGT = True;
                        if emitOnlySimpleHet:
                           #if geno[0] == "./." or geno[0] == "0/0" or geno[0] == ".":
                           #   filtGT = True;
                           anyObsAlt = False;
                           for gtt in gtTag:
                              if gtt in fmt:
                                 gtTagIdx = fmt.index(gtt);
                                 if ( geno[gtTagIdx] == "0/1" or geno[gtTagIdx] == "1/0" ):
                                    anyObsAlt = True;
                           if ( not anyObsAlt ) and emitPctRef < 0:
                              filtGT = True;
                           elif (not anyObsAlt) and random.random() > emitPctRef:
                              filtGT = True;
                        if emitOnlyNonref:
                           anyObsAlt = False;
                           for gtt in gtTag:
                              if gtt in fmt:
                                 gtTagIdx = fmt.index(gtt);
                                 genoCells = re.split("[|/]",geno[gtTagIdx])
                                 if "1" in genoCells:
                                    anyObsAlt = True;
                           if ( not anyObsAlt ) and emitPctRef < 0:
                              filtGT = True;
                           elif (not anyObsAlt) and random.random() > emitPctRef:
                              filtGT = True;
                        if emitAnyNonref:
                           anyObsAlt = False;
                           for gtt in gtTag:
                              if gtt in fmt:
                                 gtTagIdx = fmt.index(gtt);
                                 genoCells = re.split("[|/]",geno[gtTagIdx])
                                 for gg in genoCells:
                                    if gg not in [".","0"]:
                                       anyObsAlt = True;
                           if ( not anyObsAlt ) and emitPctRef < 0:
                              filtGT = True;
                           elif (not anyObsAlt) and random.random() > emitPctRef:
                              filtGT = True;
                        if emitOnlySamplesOnList:
                           if sampid not in emitOnlySamplesOnList:
                              filtGT = True;
                        if not filtGT:
                           if gtLimit > -1 and printGtCt > gtLimit:
                              filtGT = True;
                              #if printGtCt == gtLimit:
                              #   eprint("Reached limit");
                           printGtCt = printGtCt + 1;
                        if not filtGT:
                           if tableGT:
                              out.write(outline+"\t"+str(printGtCt)+"\t"+sampid);
                              for f in fmtCols:
                                 if f in fmt and len(geno) > fmt.index(f):
                                    out.write("\t"+OPENQUOTE+geno[fmt.index(f)]+QUOTECHAR);
                                 else:
                                    out.write("\t"+OPENQUOTE+NULLCHAR+QUOTECHAR);
                              if len(sampDict[sampid]) > 0:
                                 out.write("\t"+"\t".join(sampDict[sampid]));
                              out.write("\n");
                           elif flatTableGT:
                              #out.write("\t"+titleLineCells[idx]);
                              if WARN_GENO_MISSING < WARN_GENO_MISSING_LIMIT and len(geno) != len(fmt):
                                 eprint("WARNING: fmt.index(f)="+str(fmt.index(f))+", len(geno)="+str(len(geno))+"!\nFMT=\""+cells[FORMATCOLUMN]+"\"\nGENO="+genoString);
                                 WARN_GENO_MISSING = WARN_GENO_MISSING + 1;
                              if emitOnlySimpleHet or emitOnlyNonref or filterGeno:
                                 out.write("\t"+titleLineCells[idx]);
                              for f in fmtCols:
                                 if f in fmt and fmt.index(f) < len(geno):
                                    out.write("\t"+OPENQUOTE+geno[fmt.index(f)]+QUOTECHAR);
                                 #elif f in fmt:
                                 #   eprint("WARNING: fmt.index(f)="+str(fmt.index(f))+", len(geno)="+str(len(geno))+"!\nFMT=\""+cells[FORMATCOLUMN]+"\"\nGENO="+genoString);
                                 #elif f in fmt and fmt.index(f) >= len(geno) and len(geno) > 1:
                                 #   eprint("WARNING: fmt.index(f)="+str(fmt.index(f))+", len(geno)="+str(len(geno))+"!\nFMT=\""+cells[FORMATCOLUMN]+"\"\nGENO="+genoString);
                                 #   out.write("\t"+OPENQUOTE+NULLCHAR+QUOTECHAR);
                                 else:
                                    out.write("\t"+OPENQUOTE+NULLCHAR+QUOTECHAR);
                           elif emitSampleIdsOnly:
                              out.write(titleLineCells[idx]+",");
                           else :
                              out.write("\t"+titleLineCells[idx]+":"+genoString);
                     if not tableGT:
                        out.write("\n");


            else:
               out.write(outline+"\n");
            #out.write("\n");
except IOError as e:
   if e.errno == errno.EPIPE:
      exit(0);

out.close();

#WARN_GENO_MISSING = 0;
#WARN_GENO_MISSING_LIMIT = 10;

#infoMetaData[infoKey] = [0 infoKey,1 infoNum,2 infoTy,3 infoDesc,4 missingCt,5 lenMin,6 lenMax,7 min,8 max,9 examples]
#   infoWriter.write("infoTag\tnumber\ttype\tmissingCt\tlengthRange\trange\tdescription\n");
#md = [tag,num,type,desc,missinctCt,lenMin,lenMax,min,max,examples]
if writeExtendedInfoFile:
   infoWriter = open(infoFile,'w');
   infoWriter.write("infoTag\tnumber\ttype\tmissingCt\tlengthRange\trange\tdescription\n");
   for infoTag in infoCols:
      if infoTag.startswith("URL:"):
        infoWriter.write('\t'.join([infoTag,"1","URL",".",".",".","EXCEL-formatted Hyperlinked URL for "+infoTag]) + "\n");
      else:
        md = infoMetaData[infoTag]
        infoNum = md[1]
        infoTy = md[2]
        desc = md[3]
        missingCt=str(md[4])
        lenMin = md[5]
        lenMax = md[6]
        valMin = md[7]
        valMax = md[8]
        exampleSet = md[9];
        if lenMin == float('inf'):
           if infoNum == "1":
              lenRange = OPENQUOTE+"1"+QUOTECHAR;
           else:
              lenRange = OPENQUOTE+"."+QUOTECHAR;
        elif lenMin == lenMax:
           lenRange = OPENQUOTE+str(lenMin)+QUOTECHAR;
        else:
           lenRange = OPENQUOTE+str(lenMin)+"-"+str(lenMax)+QUOTECHAR;
        if valMin == float('inf'):
           valRange = OPENQUOTE+"."+QUOTECHAR;
        elif valMin == valMax:
           valRange = OPENQUOTE+str(valMin)+QUOTECHAR;
        else:
           valRange = OPENQUOTE+str(valMin)+"-"+str(valMax)+QUOTECHAR;
        if len(exampleSet) > 0:
           examples = OPENQUOTE+",".join(exampleSet)+QUOTECHAR;
        else:
           examples = valRange;
        infoWriter.write('\t'.join([infoTag,infoNum,infoTy,missingCt,lenRange,examples,desc]) + "\n");
   infoWriter.close();


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


