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

