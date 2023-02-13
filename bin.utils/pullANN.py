#!/usr/bin/env python
from __future__ import print_function

import sys,gzip,errno

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

args = sys.argv;

#+SHORTDESC="Extracts data from ANN field (from CGR VCF file)."
#+ HELPTEXT="pullANN.py utility. Extracts useful information out of the ANN field and adds it to the vcf."$'\n'
#++HELPTEXT=" "$'\n'
#+ SYNTAXTEXT="pullANN.py myVcfFile.vcf > myNewVcf.vcf"$'\n'
#++SYNTAXTEXT="pullANN.py myVcfFile.vcf.gz | gzip > myNewVcf.vcf.gz"$'\n'
#++SYNTAXTEXT="pullANN.py [options] myVcfFile.vcf.gz | gzip > myNewVcf.vcf.gz"$'\n'
#++SYNTAXTEXT="cat myVcfFile.vcf | pullANN.py [options] - > myNewVcf.vcf"$'\n'
#++SYNTAXTEXT="etc..."$'\n'
#+ PARAMS="--help: displays syntax and help info."$'\n'
#++PARAMS="--man: synonym for --help."$'\n'
#++PARAMS="--noShort: do not include short-form fields"$'\n'
#++PARAMS="--noLong: do not include long-form fields"$'\n'
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

if inf[-3:] == ".gz":
   inf = gzip.open(inf,'r');
else:
   inf = open(inf,'r');

if "--noShort" in args:
   shrtBool = False;
else:
   shrtBool = True;

if "--noLong" in args:
   longBool = False;
else:
   longBool = True;

out = sys.stdout;

lnct = 0;
isOnHeader = True;
ALEN=16;

try:
   for line in inf:
      lnct = lnct + 1;
      if line[0:2] != "##" and isOnHeader:
         if(shrtBool): out.write("##INFO=<ID=SH_HGVSc_SHORT,Number=.,Type=String,Description=\"Simple list of HGVS cDNA variant descriptions as found in the ANN field.\">\n");
         if(longBool): out.write("##INFO=<ID=SH_HGVSc_LONG,Number=.,Type=String,Description=\"List of HGVS cDNA variant descriptions as found in the ANN field. Can contain multiple entries. Each entry is labelled with the allele, the feature ID, and the HGVS cDNA variant.\">\n");
         if(shrtBool): out.write("##INFO=<ID=SH_HGVSc_SHORT_PC,Number=.,Type=String,Description=\"As SH_HGVSc_SHORT except restricted to only transcript features listed as protein coding\">\n");
         if(longBool): out.write("##INFO=<ID=SH_HGVSc_LONG_PC,Number=.,Type=String,Description=\"As SH_HGVSc_LONG except restricted to only transcript features listed as protein coding\">\n");
         if(shrtBool): out.write("##INFO=<ID=SH_HGVSp_SHORT,Number=.,Type=String,Description=\"Simple list of HGVS protein variant descriptions as found in the ANN field\">\n");
         if(longBool): out.write("##INFO=<ID=SH_HGVSp_LONG,Number=.,Type=String,Description=\"List of HGVS protein variant descriptions as found in the ANN field. Can contain multiple entries. Each entry is labelled with the allele, the feature ID, and the HGVS protein variant.\">\n");
         if(shrtBool): out.write("##INFO=<ID=SH_VARTYPE_SHORT,Number=.,Type=String,Description=\"List of variant type descriptions as found in the ANN field. Can contain multiple entries. Each entry is labelled with the allele, the feature ID, and the variant type.\">\n");
         if(longBool): out.write("##INFO=<ID=SH_VARTYPE_LONG,Number=.,Type=String,Description=\"List of variant type descriptions as found in the ANN field. Can contain multiple entries. Each entry is labelled with the allele, the feature ID, and the variant type.\">\n");
         if(longBool): out.write("##INFO=<ID=SH_VARTYPE_SHORT_PC,Number=.,Type=String,Description=\"As SH_VARTYPE_SHORT, except only for protein coding annotations/transcripts\">\n");
         if(longBool): out.write("##INFO=<ID=SH_VARTYPE_LONG_PC,Number=.,Type=String,Description=\"As SH_VARTYPE_LONG, except only for protein coding annotations/transcripts\">\n");
         isOnHeader = False;
      if line[0] == '#':
         out.write(line);
      else:
         cells = line[0:-1].split("\t");
         infocells = cells[7].split(";");
         infoPairs = [a.split('=',1) for a in infocells];
         infoNames = [a[0] for a in infoPairs];
         if not "ANN" in infoNames:
            out.write(line);
         else:
            annIdx = infoNames.index("ANN");
            annString = infoPairs[annIdx][1];
            annData = [s.strip().split("|") for s in annString.split(",")];
            cShort = [];
            cShortPC = [];
            pShort = [];
            cLong = [];
            cLongPC = [];
            pLong = [];
            tShort = [];
            tLong = [];
            tShortPC = [];
            tLongPC = [];
            for j in range(0,len(annData)):
               ad = annData[j];
               if len(ad) % ALEN != 0:
                  eprint("WARNING: bad ANN field (vcf line "+str(lnct)+"). Length is "+str(len(ad)));
               for i in range(0,len(ad)/ALEN):
                  #0-Allele, 1-Annotation, 2-Annotation_Impact, 3-Gene_Name, 4-Gene_ID, 5-Feature_Type, 6-Feature_ID, 7-Transcript_Biotype
                  #8-Rank, 9-HGVS.c, 10-HGVS.p, 11-cDNA.pos/cDNA.length, 12-CDS.pos / CDS.length, 13-AA.pos/AA.length, 14.Distance, 15.ERRORS/WARN/INFO
                  imod = i * ALEN;
                  tidx = 1 + imod;
                  cidx = 9 + imod;
                  pidx = 10 + imod;
                  feat = ad[6 + imod];
                  alle = ad[0 + imod];
                  if ad[tidx] != "":
                     tShort = tShort + [ad[tidx]];
                     tLong = tLong + [alle + "|" + feat + "|" + ad[tidx]];
                     if ad[5 + imod] == "transcript" and ad[7 + imod] == "protein_coding":
                        tShortPC = tShortPC + [ad[tidx]];
                        tLongPC = tLongPC + [alle + "|" + feat + "|" + ad[tidx]];
                  if ad[cidx] != "":
                     cLong = cLong + [alle + "|" + feat + "|" + ad[cidx]];
                     cShort = cShort + [ad[cidx]];
                     if ad[5 + imod] == "transcript" and ad[7 + imod] == "protein_coding":
                        cLongPC = cLongPC + [alle + "|" + feat + "|" + ad[cidx]];
                        cShortPC = cShortPC + [ad[cidx]];
                  if ad[pidx] != "":
                     pLong = pLong + [alle + "|" + feat + "|" + ad[pidx]];
                     pShort = pShort + [ad[pidx]];
            outInfoAdditions = [];
            if shrtBool:
               outInfoAdditions = outInfoAdditions + ["SH_HGVSc_SHORT="+(','.join(cShort) if len(cShort) > 0 else "NA")];
               outInfoAdditions = outInfoAdditions + ["SH_HGVSc_SHORT_PC="+(','.join(cShortPC) if len(cShortPC) > 0 else "NA")];
               outInfoAdditions = outInfoAdditions + ["SH_HGVSp_SHORT="+(','.join(pShort) if len(pShort) > 0 else "NA")];
               outInfoAdditions = outInfoAdditions + ["SH_VARTYPE_SHORT="+(','.join(tShort) if len(tShort) > 0 else "NA")];
               outInfoAdditions = outInfoAdditions + ["SH_VARTYPE_SHORT_PC="+(','.join(tShortPC) if len(tShortPC) > 0 else "NA")];
            if longBool:
               outInfoAdditions = outInfoAdditions + ["SH_HGVSc_LONG="+(','.join(cLong) if len(cLong) > 0 else "NA")];
               outInfoAdditions = outInfoAdditions + ["SH_HGVSc_LONG_PC="+(','.join(cLongPC) if len(cLongPC) > 0 else "NA")];
               outInfoAdditions = outInfoAdditions + ["SH_HGVSp_LONG="+(','.join(pLong) if len(pLong) > 0 else "NA")];
               outInfoAdditions = outInfoAdditions + ["SH_VARTYPE_LONG="+(','.join(tLong) if len(tLong) > 0 else "NA")];
               outInfoAdditions = outInfoAdditions + ["SH_VARTYPE_LONG_PC="+(','.join(tLongPC) if len(tLongPC) > 0 else "NA")];
            cells[7] = cells[7] + ";" + ';'.join(outInfoAdditions);
            line = '\t'.join(cells);
            out.write(line + "\n");
except IOError as e:
   if e.errno == errno.EPIPE:
      exit(0);


#
#  ../pullANN.py variants_annotated.subset.vcf.gz | less -S
#
#[C|sequence_feature|LOW|ITM2B|ENSG00000136156|topological_domain:Lumenal|ENST00000378565|protein_coding||c.565-192T>C||||||
#C|sequence_feature|LOW|ITM2B|ENSG00000136156|domain:BRICHOS|ENST00000378565|protein_coding||c.565-192T>C||||||
#C|sequence_feature|LOW|ITM2B|ENSG00000136156|disulfide_bond|ENST00000378565|protein_coding||c.565-192T>C||||||
#C|upstream_gene_variant|MODIFIER|ITM2B|ENSG00000136156|transcript|ENST00000463839|nonsense_mediated_decay||n.-2T>C|||||2532|
#C|intron_variant|MODIFIER|ITM2B|ENSG00000136156|transcript|ENST00000378565|protein_coding|4/5|c.565-192T>C||||||
#C|intron_variant|MODIFIER|ITM2B|ENSG00000136156|transcript|ENST00000378549|protein_coding|2/3|c.247-192T>C||||||
#C|intron_variant|MODIFIER|ITM2B|ENSG00000136156|transcript|ENST00000607866|nonsense_mediated_decay|4/4|n.*291-192T>C||||||]	NA
