#!/usr/bin/env Rscript 

args = commandArgs(trailingOnly=TRUE)

#args = c("/mnt/nfs/gigantor/ifs/DCEG/Projects/CCSS/steve/resources/ivs/variant_calling_intervals_120430_HG19_ExomeV3_UTR_merged_padded250bp_merged_4000parts/split/split.spans.noGap.txt",
#           "chr10","92412435")

spanFile <- args[1];
chrom <- args[2];
pos <- as.numeric( args[3] );

d <- read.table(args[1],header=F,sep='\t',stringsAsFactors=F);
names(d) <- c("partID","chrom","subpart","start","end");

which.idx <- which( d$chrom == chrom & d$start <= pos & d$end >= pos )

if(length(which.idx) == 0){
  message("ERROR: no span found!")
} else if(length(which.idx) > 1){
  message("ERROR: impossible state! Multiple spans found?")
  message( paste0(d$partID[which.idx],collapse=","))
} else {
  message(d$partID[ which.idx ])
}
