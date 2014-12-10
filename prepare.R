# Author: Matt Watts
# Date: 10 Dec 2014
# Purpose: prepare environment for RunMarxanTas R Shiny web app
#          install all the packages we need
#          preparse the spatial data structure for fast map performance

# install packages
install.packages(c("shiny","sp","maptools","PBSmapping","foreign","sqldf","vegan","labdsv","xtable","foreach","doMC"))
# you might have to install rgdal from source depending on your operating system
install.packages("rgdal",type="source")

# load packages
require(shiny)
require(sp)
require(maptools)
require(PBSmapping)
require(foreign)
require(sqldf)
require(vegan)
require(labdsv)
require(xtable)
library(foreach)
library(doMC)
library(rgdal)

# set the working directory to where your app is
setwd("/Users/matt/Documents/R/GitHub/RunMarxanTas")

# load the planning unit layer and the planning unit outline
v.pulayer <- readOGR("pulayer","pulayer")
v.puoutline <- readOGR("pulayer","tas")
# convert them
pulayer_ <<- SpatialPolygons2PolySet(v.pulayer)
puoutline <<- SpatialPolygons2PolySet(v.puoutline)
# plot them
plotPolys(pulayer_)
plotPolys(puoutline)

# load the planning unit status
pustatus <- read.csv("input/pu.dat")
pustatus_ <<- unlist(sqldf("SELECT status from pustatus"))

# save the objects to an .RData file for loading into the app at runtime
save(pulayer_,puoutline,pustatus_,
     file="pulayer/pulayer.RData")

# to launch the app in a web browser we use a command like this
shiny::runApp("/Users/matt/Documents/R/GitHub/RunMarxanTas")
