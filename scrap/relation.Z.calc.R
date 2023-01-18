








NN <- 5000000

genSamp <- function(i, N = NN){
  return(matrix(c( paste0(rep(i,N),"a"),paste0(rep(i,N),"b" ) ) ,ncol=2 ))
}

getAlle <- function(a){
   xv <- floor(runif(nrow(a),0,2)) == 1
   #out <- a;
   #out[xv,1] <- a[xv,2]
   #out[xv,2] <- a[xv,1]
   out <- a[,1];
   out[xv] <- a[xv,2];
   return(out)
}

mateSamps <- function(a,b){
   a <- getAlle(a);
   b <- getAlle(b);
   return(matrix(c(a,b),ncol=2));
}

calcZ <- function(a,b,roundTo=3,stringFmt=TRUE){
   a <- as.data.frame(a);
   b <- as.data.frame(b);
   Z0t <- a[[1]] != b[[1]] & a[[1]] != b[[2]] & a[[2]] != b[[1]] & a[[2]] != b[[2]] 
   Z2t <- ( a[[1]] == b[[1]] & a[[2]] == b[[2]] ) | 
          ( a[[1]] == b[[2]] & a[[2]] == b[[1]] )
   Z1t <- ( ! Z0t ) & ( ! Z2t )
   out <- c(mean(Z0t),mean(Z1t),mean(Z2t));
   out <- round( c( out, (out[[2]]*0.5)+out[[3]] ),roundTo);
   if(stringFmt){
     out <- sprintf(p("%.",roundTo,"f"),out)
   }
   return( out )
}

p00 <- genSamp("p00");
p01 <- genSamp("p01");
p1.sib <- mateSamps(p00,p01)


p1 <- mateSamps(p00,p01)
p2 <- genSamp("p2")
p12 <- genSamp("p12");
c1 <- mateSamps(p1,p2);
c2 <- mateSamps(p1,p2);
hs1 <- mateSamps(p1,p12);

IL1 <- genSamp("IL1");
IL2 <- genSamp("IL2");

cous1 <- mateSamps(c1,IL1);
cous2 <- mateSamps(c2,IL2);
halfcous   <- mateSamps(hs1,IL2);

IL22 <- genSamp("IL22");
IL23 <- genSamp("IL23");

sc.gp <- mateSamps(p00,p01);
sc.p  <- mateSamps(sc.gp,IL22);
sc    <- mateSamps(sc.p,IL23);

hsc1 <- mateSamps(c1,IL1);
hsc2 <- mateSamps(c2,IL1);

p3 <- genSamp("p3")
p4 <- genSamp("p4")
c3 <- mateSamps(p3,p4);
c4 <- mateSamps(p3,p4);
IL3 <- genSamp("IL3")
IL4 <- genSamp("IL4")

dc1 <- mateSamps(c1,c3);
dc2 <- mateSamps(c2,c4);

sdc1 <- mateSamps(dc1,dc2);
sdc2 <- mateSamps(dc1,dc2);

spac1 <- mateSamps(cous1,cous2)
spac2 <- mateSamps(cous1,cous2)

spas1 <- mateSamps(c1,c2);
spas2 <- mateSamps(c1,c2);

IL5 <- genSamp("IL5")
sopipos1 <- mateSamps(spas1,IL5);
sopipos2 <- mateSamps(spas1,IL5);
sopipoc1 <- mateSamps(spac1,IL5);
sopipoc2 <- mateSamps(spac1,IL5);

spahs1 <- mateSamps(hs1,c1);
spahs2 <- mateSamps(hs1,c1);

spopc1 <- mateSamps(p1,c1);
spopc2 <- mateSamps(p1,c1);

spopc3 <- mateSamps(p1,c2);
spopc4 <- mateSamps(p1,c2);


RL <- 
list(
c("parent-offspring",calcZ(p1,c1)),
c("siblings",calcZ(c1,c2)),
c("cousins",calcZ(cous1,cous2)),
c("aunt-niece",calcZ(cous1,c2)),
c("halfsiblings",calcZ(c1,hs1)),
c("halfsiblings/cousins",calcZ(hsc1,hsc2)),
c("siblings, parents are cousins",calcZ(spac1,spac2)),
c("double cousins",calcZ(dc1,dc2)),
c("parent-offspring, parents are cousins",calcZ(cous1,spac1)),
c("siblings, parents are double cousins",calcZ(sdc1,sdc2)),
c("parent-offspring, parents are double cousins",calcZ(dc1,sdc1)),

c("g.parent-g.child",calcZ(p1,cous1)),
c("g.parent-g.child, parents are cousins",calcZ(p1,spac1)),

c("second cousins",calcZ(sc,cous1)),
c("half cousins",calcZ(halfcous,cous1)),

c("gg.parent-gg.child",calcZ(p00,cous1)),
c("g.aunt-g.niece",calcZ(p1.sib,cous1)),

c("siblings, parents are half-siblings",calcZ(spahs1 ,spahs2)),
c("parent-offspring, parents are half-siblings",calcZ(c1 ,spahs2)),

c("siblings, parents are siblings",calcZ(spas1,spas2)),
c("parent-offspring, parents are siblings",calcZ(c1,spas1)),
c("siblings, one parent is child of cousins",calcZ(sopipos1 ,sopipos2 )),
c("siblings, one parent is child of siblings",calcZ(sopipoc1 ,sopipoc2 )),

c("parent-child, parent is child of cousins",calcZ(spas1 ,sopipos2 )),
c("parent-child, parent is child of siblings",calcZ(spac1 ,sopipoc1 )),

c("siblings, parent is parent-child",calcZ(spopc1 ,spopc2 )),
c("parent-child, parent is parent-child, gparent",calcZ(p1 ,spopc1 )),
c("parent-child, parent is parent-child, child",  calcZ(c1 ,spopc1 ))
)

rd <- do.call(rbind.data.frame,RL)
names(rd) <- c("relation","Z0","Z1","Z2","PI_HAT")
rd

write.table(rd,file="Z.distribution.by.pedigree.v010.txt",row.names=F,col.names=T,sep='\t',quote=F);



c("grandniece-greataunt",calcZ(p1,sc)),

####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################


#parent-offspring:
calcZ(p1,c1)
#siblings:
calcZ(c1,c2)
#cousins:
calcZ(cous1,cous2)
#uncle-nephew
calcZ(cous1,c2)
calcZ(cous2,c1)
#half sibling:
calcZ(c1,hs1)
#half-sibling / cousin:
calcZ(hsc1,hsc2)
#siblings, parents are cousins:
calcZ(spac1,spac2)
#siblings, parents are siblings:
calcZ(spas1,spas2)
#double-cousins (cousins on both sides)
calcZ(dc1,dc2)
#parent-offspring where parents are siblings:
calcZ(c1,spas1)
#parent-offspring where parents are cousins:
calcZ(cous1,spac1)
#





