#!/usr/bin/env Rscript

#TODO-WRITE-HELPDOC

args = commandArgs(trailingOnly=TRUE)

infile <- args[length(args)-1];
outfile <- args[length(args)];

d <- as.matrix( read.table(infile,header=F,sep='\t',stringsAsFactors=F) );

write.table(t(d),file=outfile,row.names=F,col.names=F,quote=F,sep='\t');

