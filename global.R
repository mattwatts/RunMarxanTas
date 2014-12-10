# GitHub RunMarxanTas

library(shiny)
library(PBSmapping)
library(maptools)
library(sp)

cat(paste0("hello\n"))
cat(paste0(getwd(),"\n"))

sMarxanDir <- getwd()

inputdat <- readLines(paste(sMarxanDir,"/input.dat",sep=""))

iParam <- which(regexpr("NUMREPS",inputdat)==1)
iNUMREPS <<- as.integer(unlist(strsplit(inputdat[iParam], split=" "))[2])

irefreshinput <<- 0
isavetargetspf <<- 0
