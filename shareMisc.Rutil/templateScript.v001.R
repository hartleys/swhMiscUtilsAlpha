#I always set this option:
options(stringsAsFactors=FALSE);

#This line autodetects if you are a windows machine with T drive access
#    or on CCAD or on Helix/Biowulf and points to the right file location:
SHAREMISC.DIR <- ifelse(dir.exists("T:/"),"T:/shared/hartleys/software/shareMisc/",
                 ifelse(dir.exists("/mnt/nfs/gigantor/ifs/Shared/hartleys/software/shareMisc/"),"/mnt/nfs/gigantor/ifs/Shared/hartleys/software/shareMisc/",
                 ifelse(dir.exists("/data/hartleys/pub/software/shareMisc/"),"/data/hartleys/pub/software/shareMisc/",stop("ERROR: STEVE HARTLEY SCRIPTS NOT FOUND!"))));
#Load all my functions:
source(paste0( SHAREMISC.DIR,"/shareMisc.Rutil/misc.R"));

TDRIVE <- ""
if(file.exists("T:\\")){
  TDRIVE <- "T:\\"
}



###########################################################################################################################
###########################################################################################################################
###########################################################################################################################
###########################################################################################################################
###########################################################################################################################
###########################################################################################################################
---------------------------------------------------------------------------------------------------------------------------
---------------------------------------------------------------------------------------------------------------------------
---------------------------------------------------------------------------------------------------------------------------





















