## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
library(knitr)
library(PersonAlyticsPower)
papv <- packageVersion("PersonAlyticsPower")


depends <- function()
{
  d <- scan('../DESCRIPTION', what = 'character', sep='\n')
  wd <- which( grepl("*Depends*", d) )
  wr <- which( grepl("RoxygenNote*", d) )
  pkgs <- c(d[(wd+1):(wr-1)])
  pkgs <- gsub( " ", "", pkgs)
  pkgs <- gsub( "\t", "", pkgs)
  pkgs <- gsub( ",", "", pkgs)
  pkgs
}
pkgs <- depends()
pkgss <-  matrix( unlist( strsplit(pkgs, "\\(") ), ncol=2, byrow=TRUE)[,1]
pkgss <- pkgss[pkgss!='PersonAlytics']

## ----eval=FALSE----------------------------------------------------------
#  install.deps <- function(pkg)
#  {
#    new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
#    old.pkg <- pkg[ (pkg %in% installed.packages()[, "Package"])]
#    if (length(new.pkg)) install.packages(new.pkg, dependencies = TRUE,
#                                          repos = "https://cran.rstudio.com/")
#    if (length(old.pkg)) update.packages(old.pkg, dependencies = TRUE)
#  }

## ---- echo=FALSE, comment=NA---------------------------------------------
cat(paste("install.deps(c(", paste("'", pkgss[-1], "'", collapse=', ', sep=""), "))", 
          sep=""), "\n")

## ---- message=FALSE------------------------------------------------------
library(PersonAlyticsPower)

## ---- eval=FALSE---------------------------------------------------------
#  ?ICTpower

## ------------------------------------------------------------------------
myPolyICT <- polyICT$new(
  groups            = c(group1=10, group2=10)                 ,
  phases            = makePhase(nObsPerPhase = c(10, 20, 10)) ,
  propErrVar        = c(randFx=.5, res=.25, mserr=.25)        ,
  randFxOrder       = 1                                       ,
  randFxCor         = 0.2                                     ,
  randFxVar         = c(1, 1)                                 ,
  error             = armaErr$new()                           ,
  merror            = armaErr$new(list())                     ,
  ySD               = 15                                      ,
  yMean             = 100                                     ,
  )

## ------------------------------------------------------------------------
myPolyICT

## ------------------------------------------------------------------------
myPolyICT$inputMat

## ----eval=FALSE----------------------------------------------------------
#  edit(myPolyICT$inputMat)

## ----echo=FALSE----------------------------------------------------------
myPolyICT$inputMat[myPolyICT$inputMat$Phase=='phase2' &
  myPolyICT$inputMat$Group=='group1', 'Mean0'] <- .3
myPolyICT$inputMat[myPolyICT$inputMat$Phase=='phase3' &
  myPolyICT$inputMat$Group=='group1', 'Mean0'] <- .3
myPolyICT$inputMat[myPolyICT$inputMat$Phase=='phase3' &
  myPolyICT$inputMat$Group=='group1', 'Mean1'] <- -.6

## ------------------------------------------------------------------------
myPolyICT$inputMat

## ---- fig.width = 6, fig.height=5----------------------------------------
myPolyICT$designCheck(ylim=c(75,125))

## ------------------------------------------------------------------------
ICTpower(outFile = c('twoGroupThreePhase', 'csv') ,
         design  = myPolyICT                      , 
         B       = 10                             ,
         prompt  = FALSE                          ) 

## ------------------------------------------------------------------------

# use a deep clone to prevent `myPolyICTnonPar` from updating `myPolyICT`
myPolyICTnonPar <- myPolyICT$clone(deep=TRUE)

# update the groups to be 500 each
myPolyICTnonPar$update(groups=c(group1=500, group2=500))

# simulate and save the data
Data <- myPolyICTnonPar$makeData()
save(Data, file = "Data.RData")

## ------------------------------------------------------------------------
ICTpower(outFile     = c("npbsTest", "csv") ,
         B           = 10                   ,
         dataFile    = "Data.RData"         ,
         sampleSizes = c(25,25)             ,
         prompt      = FALSE                )

## ------------------------------------------------------------------------
names(myPolyICT)

## ------------------------------------------------------------------------
myPolyICT$randFxFamParms
myPolyICT$randFxFamParms <- list(mu=.3, sigma=.8)
myPolyICT$randFxFamParms

## ------------------------------------------------------------------------
myPolyICT$randFxCorMat
myPolyICT$randFxCorMat$phase1$group1[2,1] <- 
  myPolyICT$randFxCorMat$phase1$group1[1,2] <- .45
myPolyICT$randFxCorMat

## ------------------------------------------------------------------------
# groups
myPolyICT$update(groups=c(group1=120, group2=15))
myPolyICT$groups

# phases
myPolyICT$update(phases=makePhase(c(10,10,10)))
myPolyICT$phases

# propErrVar - not reccomended for $update, use edit(myPolyICT$inputMat) instead

# randFxOrder - not reccomended for $update, create a new polyICT object instead

# randFxCor - not reccomended for $update, edit myPolyICT$randFxCorMat instead

# randFxVar - not reccomended for $update, use edit(myPolyICT$inputMat) instead

# error - not reccomended for $update, create a new polyICT object instead

# merror - not reccomended for $update, create a new polyICT object instead

# ySD
myPolyICT$update(ySD = 1, yMean = 0)
myPolyICT$ySD
myPolyICT$yMean


## ------------------------------------------------------------------------
ICTpower(outFile     = c("npbsFPCtest", "csv") ,
         B           = 10                      ,
         dataFile    = "Data.RData"            ,
         sampleSizes = c(25,25)                ,
         fpc         = length(table(Data$id))  ,
         prompt      = FALSE                   )

## ---- echo=FALSE, message=FALSE------------------------------------------
# clean up
txts <- dir(getwd(), glob2rx("*.txt"))
csvs <- dir(getwd(), glob2rx("*.csv"))
file.remove("Data.RData", txts, csvs)

