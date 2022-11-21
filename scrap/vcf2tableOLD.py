#!/usr/bin/env python
from __future__ import print_function

import sys,gzip,errno

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

args = sys.argv;

#+SHORTDESC="Converts VCF to tab-delimited table"
#+ HELPTEXT="vcf2table.py utility. Takes the info data from a VCF and makes a table. For testing purposes."$'\n'
#++HELPTEXT=" "$'\n'
#+ SYNTAXTEXT="vcf2table.py [options] myVcf.vcf.gz"$'\n'
#+ PARAMS="--help: displays syntax and help info."$'\n'
#++PARAMS="--man: synonym for --help."$'\n'
#++PARAMS="--infoColumnList <filename>: a comma-delimited list of INFO columns you want extracted, in order."$'\n'
#++PARAMS="--keepGT: Flag, keep genotype columns."$'\n'
#++PARAMS="--vcfColumnList <filename>: a comma-delimited list of INFO columns you want extracted, in order."$'\n'
#++PARAMS="--infoColumnList <list>: a comma-delimited list of primary VCF columns you want included, in order. By default all are included."$'\n'
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


if "--infoColumnList" in args:
   idx = args.index("--infoColumnList");
   infoCols = args.pop(idx+1).split(",");
   args.pop(idx);
   presetInfoCols = True;

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

infoFlagSet = set();

try:
   for line in inf:
      lnct = lnct + 1;
      if lnct % 1000 == 0 and not isOnHeader:
         eprint(".");
      
      if len(line) == 1:
         out.write("\n");
      elif line[0:2] != "##" and isOnHeader:
         out.write("CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\t"+'\t'.join(infoCols)+"\n");
         isOnHeader = False;
         eprint("Finished with header...")
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
      elif line[0] != "#":
         cells = line[0:-1].split("\t");
         infocells = cells[7].split(";");
         infoPairs = [a.split('=',1) for a in infocells];
         infoNames = [a[0] for a in infoPairs];
         out.write('\t'.join(cells[0:7]));
         
         for c in infoCols:
            if c in infoNames:
               if c in infoFlagSet:
                  out.write("\t"+QUOTECHAR+"TRUE"+QUOTECHAR);
               else:
                  idx = infoNames.index(c)
                  if len(infoPairs[idx]) == 1:
                     out.write("\t"+QUOTECHAR+BLANKCHAR+QUOTECHAR);
                  else:
                     out.write("\t"+QUOTECHAR+infoPairs[idx][1]+QUOTECHAR);
            else:
               if c in infoFlagSet:
                  out.write("\t"+QUOTECHAR+"FALSE"+QUOTECHAR);
               else:
                  out.write("\t"+QUOTECHAR+NULLCHAR+QUOTECHAR);
         if(keepGT):
            out.write( "\t"+ "\t".join(cells[8:]) );
         out.write("\n");
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

