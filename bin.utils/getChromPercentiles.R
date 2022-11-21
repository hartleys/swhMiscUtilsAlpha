#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

chromList = args[length(args)];

d <- read.table(chromList,header=F,sep='\t',stringsAsFactors=F);
names(d) <- c("chrom","len");
d$MB <- d$len / 1000000;
d$sumMB <- cumsum( d$MB );
totMB <- sum(d$MB)
d$pct  <- d$sumMB / totMB

pct.breaks <- 0:39 / 40;
break.len <- 1 / 40;

pct.breaks <- lapply( pct.breaks, function(pct){
  mb <- pct * totMB;
  if( (d$sumMB > mb)[[1]] ){
    chrom.idx <- 1;
    mb.rem <- mb;
  } else {
    chrom.idx <- which( d$sumMB > mb )[[1]];
    mb.rem <- mb - d$sumMB[[chrom.idx - 1]];
  }
  bp.rem <-mb.rem * 1000000
  bp.rem.string <- format( mb.rem * 1000000, big.mark = ",",scientific=F)
  out <- list(pct = pct,chrom = d$chrom[chrom.idx],bp = bp.rem,bp.fmt = bp.rem.string, chrom.len = format(d$len[[chrom.idx]], big.mark = ",",scientific=F));
  #print(data.frame(field=names(out),V = unlist(out)))
  out;
})
pct.breaks <- do.call(rbind.data.frame,pct.breaks)


pct.breaks

chr5:153,619,217


