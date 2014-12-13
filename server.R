# Author: Matt Watts
# Date: 10 Dec 2014
# Purpose: RunMarxanTas web app server.R

require(shiny)
require(sp)
require(maptools)
require(PBSmapping)
require(foreign)
require(sqldf)
require(vegan)
require(labdsv)
require(xtable)
require(lpSolveAPI)
library(foreach)
library(doMC)

registerDoMC(10)  # the number of CPU cores 

cat(paste0(sMarxanDir,"\n"))

# initialise objects from R binary files
load(file=paste0(sMarxanDir,"/pulayer/pulayer.Rdata"))
puoutline <<- puoutline
pulayer_ <<- pulayer_
pustatus_ <<- pustatus_

PrepareDisplay <- function()
{
    cat("PrepareDisplay\n")

    # prepare the map: pulayer object
    pulayer <<- pulayer_
    pu_table <<- read.dbf(paste(sMarxanDir,"/pulayer/pulayer.dbf",sep=""))
    
    # prepare the planning unit status object
    pustatus <<- pustatus_

    # prepare the cluster analysis objects and map objects for ILP & LP
    solutions_raw<-read.table(paste0(sMarxanDir,"/output/output_solutionsmatrix.csv"),header=TRUE, row.name=1, sep=",")
    # move best solution to start of table
    thetable <- read.csv(paste0(sMarxanDir,"/output/output_sum.csv"))
    thetable <- round(sqldf("SELECT Score, Cost, Planning_Units, Penalty, Shortfall, Missing_Values, MPM from thetable"))
    iBest <- which.min(thetable[,1])
    Best <- solutions_raw[iBest,]
    solutions_raw <- solutions_raw[-iBest,]
    solutions_join <- rbind(Best,solutions_raw)

    rownames(solutions_join) <- c(paste0("S",iBest," (Best)"),row.names(solutions_raw))
    plotlabels <<- c(paste0("S",iBest," (Best)"),row.names(solutions_raw))
    
    solutions <- unique(solutions_join)
    iUniqueSolutions <- dim(solutions)[1]
    nmdscolours <- rep("black",each = iUniqueSolutions)  
    nmdscolours[1] <- "blue"
    nmdscolours <<- nmdscolours
    soldist<-vegdist(solutions,distance="bray")
    sol.mds<<-nmds(soldist,2)
    h<<-hclust(soldist, method="complete")
}

PrepareDisplay()

GetOutputFileext <- function(sMarxanDir,sParam)
# For the specified Marxan output file, return the file extension (.csv or .txt)
# Scan input.dat for the parameter,
# if value = 1, .dat, tab delimited, no header
# if value = 2, .txt, comma delimited (Qmarxan sets this value)
# if value = 3, .csv, comma delimited
{
  inputdat <- readLines(paste(sMarxanDir,"/input.dat",sep=""))
  iParam <- which(regexpr(sParam,inputdat)==1)
  
  iValue <- as.integer(unlist(strsplit(inputdat[iParam], split=" "))[2])
  
  if (iValue == 1)
  {
    return(".dat")
  }
  if (iValue == 2)
  {
    return(".txt")
  }
  if (iValue == 3)
  {
    return(".csv")
  }
}

GenerateSolnFilename <- function(iRunNumber,sMarxanDir)
{
  sFilename <- paste(sMarxanDir,"/output/output_r",sep="")  
  iPadding <- 5 - nchar(as.character(iRunNumber))
  if (iPadding > 0)
  {
    for (i in 1:iPadding)
    {
      sFilename <- paste(sFilename,"0",sep="")
    }
  }
  sFilename <- paste(sFilename,iRunNumber,GetOutputFileext(sMarxanDir,"SAVERUN"),sep="")  
}

ImportOutputsCsvToShpDbf <- function(sPuShapeFileDbf, sMarxanDir, iNumberOfRuns, sPUID)
# Imports the relevant contents of output files to the planning unit shape file dbf.
{
  # load and prepare pu_table
  pu_table <- read.dbf(sPuShapeFileDbf)
  pu_table <- sqldf(paste("SELECT ", sPUID, " from pu_table",sep=""))
  colnames(pu_table)[1] <- "PUID"
                    
  pu_table$PUID <- as.integer(pu_table$PUID)
  
  # load and prepare ssoln_table
  ssoln_table <- read.csv(paste(sMarxanDir,"/output/output_ssoln",GetOutputFileext(sMarxanDir,"SAVESUMSOLN"),sep=""))
  colnames(ssoln_table)[1] <- "PUID"
  colnames(ssoln_table)[2] <- "SSOLN2"
  ssoln_table$SSOLN1 <- as.integer(iNumberOfRuns - ssoln_table$SSOLN2)
  ssoln_table$SSOLN2 <- as.integer(ssoln_table$SSOLN2)
  
  # join pu_table and ssoln_table
  pu_table <- sqldf("SELECT * from pu_table LEFT JOIN ssoln_table USING(PUID)")
  
  # load and prepare best_table
  best_table <- read.csv(paste(sMarxanDir,"/output/output_best",GetOutputFileext(sMarxanDir,"SAVEBEST"),sep=""))
  best_table$BESTSOLN <- as.integer(best_table$SOLUTION + 1)
  best_table <- sqldf("SELECT PUID, BESTSOLN from best_table")
  
  # join pu_table and best_table
  pu_table <- sqldf("SELECT * from pu_table LEFT JOIN best_table USING(PUID)")
  
  # save the new pu_table
  colnames(pu_table)[1] <- sPUID
  write.dbf(pu_table,sPuShapeFileDbf)  
}

labelCol <- function(x)
{
  thetable <- read.csv(paste0(sMarxanDir,"/output/output_sum.csv"))
  thetable <- round(sqldf("SELECT Score from thetable"))
  iBest <- which.min(thetable[,1])
  if (is.leaf(x))
  {
    ## fetch label
    a <- attributes(x)
    label <- attr(x, "label") 
    colour <- "black"
    if (label == paste0("S",iBest," (Best)")) { colour <- "blue" } #"green" }
    attr(x, "nodePar") <- c(a$nodePar, lab.col = colour)
  }
  return(x)
}

PadInt <- function(iRunNumber)
{
  sFilename <- ""
  iPadding <- 5 - nchar(as.character(iRunNumber))
  if (iPadding > 0)
  {
    for (i in 1:iPadding)
    {
      sFilename <- paste0(sFilename,"0")
    }
  }
  sFilename <- paste0(sFilename,iRunNumber)  
}

JoinParallelResults <- function()
{
    # combine the summary tables
    sumtable <- c()
    for (i in 1:10)
    {
        sumtable_ <- read.csv(paste0(sMarxanDir,"/output/output",i,"_sum.csv"))
        sumtable <- rbind(sumtable,sumtable_)
    }
    for (i in 1:100)
    {
        sumtable[i,1] <- i
    }
    write.csv(sumtable,
              paste0(sMarxanDir,"/output/output_sum.csv"),
              quote=FALSE,row.names=FALSE)

    # detect best solution
    iBest <- which(sumtable[,2]==min(sumtable[,2]))
    if (length(iBest) > 0)
    {
        iBest <- iBest[1]
    }
    
    # rename mv files and solution files
    iSol <- 0
    for (i in 1:10)
    {
        for (j in 1:10)
        {
            iSol <- iSol + 1
        
            file.rename(paste0(sMarxanDir,"/output/output",i,"_mv",PadInt(j),".csv"),
                        paste0(sMarxanDir,"/output/output_mv",PadInt(iSol),".csv"))
        
            file.rename(paste0(sMarxanDir,"/output/output",i,"_r",PadInt(j),".csv"),
                        paste0(sMarxanDir,"/output/output_r",PadInt(iSol),".csv"))
        }
    }

    # copy _mvbest and _best files
    file.copy(paste0(sMarxanDir,"/output/output_mv",PadInt(iBest),".csv"),
              paste0(sMarxanDir,"/output/output_mvbest.csv"),
              overwrite=TRUE)
    file.copy(paste0(sMarxanDir,"/output/output_r",PadInt(iBest),".csv"),
              paste0(sMarxanDir,"/output/output_best.csv"),
              overwrite=TRUE)
    
    # join ssoln files
    ssolntable <- read.csv(paste0(sMarxanDir,"/output/output",i,"_ssoln.csv"))
    colnames(ssolntable)[2] <- "SS1"
    for (i in 2:10)
    {
        ssolntable_ <- read.csv(paste0(sMarxanDir,"/output/output",i,"_ssoln.csv"))
        ssolntable <- sqldf("SELECT * from ssolntable LEFT JOIN ssolntable_ USING(planning_unit)")
        colnames(ssolntable)[ncol(ssolntable)] <- paste0("SS",i)
    }
    ssolntable$number <- ssolntable$SS1 + ssolntable$SS2 + ssolntable$SS3 + ssolntable$SS4 + ssolntable$SS5 + ssolntable$SS6 + ssolntable$SS7 + ssolntable$SS8 + ssolntable$SS9 + ssolntable$SS10
    ssolntable <- sqldf("SELECT planning_unit, number from ssolntable")
    write.csv(ssolntable,
                  paste0(sMarxanDir,"/output/output_ssoln.csv"),
                  quote=FALSE,row.names=FALSE)

    # join cluster files: text parse
    outfile <- file(paste0(sMarxanDir,"/output/output_solutionsmatrix.csv"),"w")
    iRow <- 0
    for (i in 1:10)
    {
        infile <- file(paste0(sMarxanDir,"/output/output",i,"_solutionsmatrix.csv"),"r")
        # read header row
        sLine <- readLines(con=infile,n=1)
  
        # write header row if i == 1
        if (i == 1)
        {
            write(sLine,file=outfile)
        }
    
        for (j in 1:10)
        {
            sLine <- readLines(con=infile,n=1)
            iLen <- nchar(sLine)
            if (j < 10)
            {
                # S1..S9 : remove 3 chars
                sLine <- substr(sLine, 4, iLen)          
            } else {
                # S10 : remove 4 chars
                sLine <- substr(sLine, 5, iLen)
            }
            iRow <- iRow + 1
            write(paste0("S",iRow,",",sLine),file=outfile,append=TRUE)
        }
        close(infile)
    }
    close(outfile)
}

shinyServer(function(input, output, session) {

    observe({
        rprop <<- input$prop
    })

    observe({
        rspf <<- input$spf
    })
    
    observe({
        sfeature <<- input$feature
        
        # change the target text control for the selected feature
        specdat <- read.csv(paste(sMarxanDir,"/input/spec.dat",sep=""),stringsAsFactors=FALSE)
        for (j in 1:nrow(specdat))
        {
            if (specdat[j,4] == sfeature)
            {
                updateNumericInput(session, "prop", value = specdat[j,2])
                updateNumericInput(session, "spf", value = specdat[j,3])
            }
        }
    })

    observe({
        input$savetargetspf
        cat("savetargetspf\n")
        
        rtarget <- rprop
        if (rtarget < 0)
             { rtarget <- 0 }
        if (rtarget > 1)
            { rtarget <- 1 }
        if (rspf < 0)
            { rspf <- 0 }

        # save target/spf to spec.dat
        specdat <- read.csv(paste(sMarxanDir,"/input/spec.dat",sep=""),stringsAsFactors=FALSE)
        # change the value only for the row with name == input$feature
        for (j in 1:nrow(specdat))
        {
            if (sfeature == "All features")
            {
                cat("saving all features...\n")
                
                specdat[j,2] <- rtarget
                specdat[j,3] <- rspf
                
            } else {
                if (specdat[j,4] == sfeature)
                {
                    specdat[j,2] <- rtarget
                    specdat[j,3] <- rspf
                }
            }
        }
        write.csv(specdat,paste0(sMarxanDir,"/input/spec.dat"),quote=FALSE,row.names=FALSE)
    })

    runmarxan <- reactive({
        if (input$mrun == 0)
        {
            imrun <<- 0
            cat("init mrun\n")
        }
        else
        {
            if (input$mrun > imrun)
            {
                imrun <<- input$mrun
                cat("mrun incremented\n")

                # BLM parameter
                inputdat <- readLines(paste0(sMarxanDir,"/input.dat"))
                iBLMparam <- which(regexpr("BLM",inputdat)==1)
                inputdat[iBLMparam] <- paste0("BLM ",input$blm)
                
                randomseeds <- round(runif(10)*100000)
                
                # run Marxan
                foreach(i=1:10) %dopar%
                {
                    # set parameters for multi core
                    iINPUTDIRparam <- which(regexpr("INPUTDIR",inputdat)==1)
                    iOUTPUTDIRparam <- which(regexpr("OUTPUTDIR",inputdat)==1)
                    iSCENNAMEparam <- which(regexpr("SCENNAME",inputdat)==1)
                    iNUMREPSparam <- which(regexpr("NUMREPS",inputdat)==1)
                    iRANDSEEDparam <- which(regexpr("RANDSEED",inputdat)==1)
                    inputdat[iINPUTDIRparam] <- paste0("INPUTDIR ",sMarxanDir,"/input")
                    inputdat[iOUTPUTDIRparam] <- paste0("OUTPUTDIR ",sMarxanDir,"/output")
                    inputdat[iSCENNAMEparam] <- paste0("SCENNAME output",i)         
                    inputdat[iNUMREPSparam] <- "NUMREPS 10"
                    inputdat[iRANDSEEDparam] <- paste0("RANDSEED ",randomseeds[i])
                    writeLines(inputdat,paste0(sMarxanDir,"/core",i,"/input.dat"))
                
                    cat(paste0("getwd ",getwd(),"\n"))
                    setwd(paste0(sMarxanDir,"/core",i))
                    cat(paste0("getwd ",getwd(),"\n"))
                    cat(paste0(".Platform$pkgType ",.Platform$pkgType,"\n"))
                    if ((.Platform$pkgType == "mac.binary") || (.Platform$pkgType == "mac.binary.mavericks"))
                    {
                        cat("mac\n")
                        system("./MarOpt_v243_Mac64 -s")
                    } else {
                        cat("linux\n")
                        system("./MarOpt_v243_Linux64 -s")
                    }
                }
                
                JoinParallelResults()
                
                # write results to pu dbf
                ImportOutputsCsvToShpDbf(paste0(sMarxanDir,"/pulayer/pulayer.dbf"),
                                         sMarxanDir,iNUMREPS,"PU_ID")
                
                # fetch the results
                PrepareDisplay()

                # trigger a refresh of the UI
                irefreshinput <<- irefreshinput + 1
                updateNumericInput(session, "refreshinput", value = irefreshinput)
            }
        }

        return(as.character(input$mrun))
    })

    outputmap <- reactive({

        input$refreshinput
        
        if (input$map == "ilpmap")
        {
            # like bestmap
            greenramp <- colorRampPalette(c("white","green"))(2)
            colours <- rep(greenramp[1],length(ILPvalues))
            for (j in 1:length(ILPvalues))
            {
                if (pustatus[j] == 2)
                {
                    colours[j] <- "#40E0D0" # Turquoise
                } else {
                    if (pustatus[j] == 3)
                    {
                        colours[j] <- "grey"
                    } else {
                        colours[j] <- greenramp[ILPvalues[j]]
                   }
                }
            }
            plotPolys(pulayer,col=colours,axes=FALSE,border=NA,cex.lab=0.1,cex.axis=0.1)
        }
        if (input$map == "ssolnNmap")
        {
            values <- sqldf(paste("SELECT SSOLN2 from pu_table",sep=""))
            blueramp <- colorRampPalette(c("white","blue"))(16)
            colours <- rep(blueramp[1],nrow(values))
            for (j in 1:nrow(values))
            {
                if (pustatus[j] == 2)
                {
                    colours[j] <- "#40E0D0" # Turquoise
                } else {
                    if (pustatus[j] == 3)
                    {
                        colours[j] <- "grey"
                    } else {
                        colours[j] <- blueramp[round(15 / 100 * values[j,])+1]
                    }
                }
            }
            plotPolys(pulayer,col=colours,axes=FALSE,border=NA,cex.lab=0.1,cex.axis=0.1)
        }

        if (input$map == "bestmap")
        {
            values <- sqldf("SELECT BESTSOLN from pu_table")
            greenramp <- colorRampPalette(c("white","green"))(2)
            colours <- rep(greenramp[1],nrow(values))
            for (j in 1:nrow(values))
            {
                if (pustatus[j] == 2)
                {
                    colours[j] <- "#40E0D0" # Turquoise
                } else {
                    if (pustatus[j] == 3)
                    {
                        colours[j] <- "grey"
                    } else {
                        colours[j] <- greenramp[values[j,]]
                    }
                }
            }
            plotPolys(pulayer,col=colours,axes=FALSE,border=NA,cex.lab=0.1,cex.axis=0.1)
        }

        if (input$map == "runMmap")
        {
            solnX_table <- read.csv(GenerateSolnFilename(input$m,sMarxanDir))
            values <- sqldf(paste("SELECT SOLUTION from solnX_table",sep="")) + 1
            greenramp <- colorRampPalette(c("white","green"))(2)
            colours <- rep(greenramp[1],nrow(values))
            for (j in 1:nrow(values))
            {
                if (pustatus[j] == 2)
                {
                    colours[j] <- "#40E0D0" # Turquoise
                } else {
                    if (pustatus[j] == 3)
                    {
                        colours[j] <- "grey"
                    } else {
                        colours[j] <- greenramp[values[j,]]
                    }
                }
            }
            plotPolys(pulayer,col=colours,axes=FALSE,border=NA,cex.lab=0.1,cex.axis=0.1)
        }

        addLines(puoutline,col="black")		
    })

    outputtable <- reactive({

        input$refreshinput
        input$savetargetspf
        
        if (input$table == "spec")
        {
            thetable <- read.csv(paste0(sMarxanDir,"/input/spec.dat"))
        }
        if (input$table == "sumtable")
        {
            sFilename <- paste(sMarxanDir,"/output/output_sum.csv",sep="")
            thetable <- read.csv(sFilename)
            thetable <- round(sqldf("SELECT Score, Cost, Planning_Units, Penalty, Shortfall from thetable"))
            iBest <- which.min(thetable[,1])
            Run <- c()
            for (j in 1:nrow(thetable))
            {
                if (j == iBest)
                {
                    Run <- c(Run,"Best")
                } else {
                    Run <- c(Run,j)
                }
            }

            thetable <- cbind(Run,thetable)#,
            thetable$Run <- as.character(thetable$Run)
            thetable$Run <- as.character(thetable$Run)
            for (j in 1:6)
            {
                thetable[iBest,j] <- HTML(paste0("<FONT COLOR='blue'>",thetable[iBest,j],"</FONT>"))
            }
        }
        if (input$table == "mvbesttable")
        {
                sFilename <- paste(sMarxanDir,"/output/output_sum.csv",sep="")
                thetable <- read.csv(sFilename)
                thetable <- round(sqldf("SELECT Score, Cost, Planning_Units, Penalty, Shortfall from thetable"))
                iBest <- which.min(thetable[,1])

                sFilename <- paste(sMarxanDir,"/output/output_mv",sep="")
                iPadding <- 5 - nchar(as.character(iBest))
                if (iPadding > 0)
                {
                        for (i in 1:iPadding)
                        {
                                sFilename <- paste(sFilename,"0",sep="")
                        }
                }
                sFilename <- paste(sFilename,iBest,".csv",sep="")

                thetable <- read.csv(sFilename,stringsAsFactors=FALSE)
                # sort the table the way spec.dat is ordered
                tableorder <- c(25,24,23,22,21,20,19,18,17,16,15,14,13,12,11,10,9,8,7,6,5,4,3,2,1)
                thetable <- thetable[tableorder,]
                # select just fields we want
                colnames(thetable)[4] <- "AmountHeld"
                colnames(thetable)[9] <- "TargetMet"
                thetable <- sqldf(paste0("SELECT Target, AmountHeld, TargetMet from thetable"))
                # load feature name from spec.dat
                specdat <- read.csv(paste0(sMarxanDir,"/input/spec.dat"),stringsAsFactors=FALSE)
                name <- sqldf(paste0("SELECT name from specdat"))
                # join the name
                thetable <- cbind(name,thetable)
                # select the fields we want from total areas file
                tafile <- read.csv(paste0(sMarxanDir,"/core1/MarOptTotalAreas.csv"))
                tafile <- tafile[tableorder,]
                tafile <- sqldf(paste0("SELECT totalarea,reservedarea from tafile"))
                # compute target gap
                targetgap <- rep(0,each=nrow(thetable))
                for (i in 1:nrow(thetable))
                {
                    if (thetable$Target[i] > 0)
                    {
                        if (thetable$AmountHeld[i] < thetable$Target[i])
                        {
                            targetgap[i] <- thetable$Target[i] - thetable$AmountHeld[i]
                        }
                    }
                }
                # join and tidy the table
                thetable <- cbind(thetable,tafile,targetgap)
                thetable <- sqldf(paste0("SELECT name, totalarea, reservedarea, Target, AmountHeld, TargetMet, targetgap from thetable"))
                colnames(thetable)[2] <- "Total"
                colnames(thetable)[3] <- "Reserved"
                colnames(thetable)[7] <- "TargetGap"
                # colour code features that have not met targets
                for (i in 1:nrow(thetable))
                {
                    if (thetable[i,6] == "no")
                    {
                        for (j in 1:7)
                        {
                            thetable[i,j] <- HTML(paste0("<FONT COLOR='blue'>",thetable[i,j],"</FONT>"))
                        }
                    } else {
                        for (j in 1:7)
                        {
                            thetable[i,j] <- HTML(paste0("<FONT COLOR='black'>",thetable[i,j],"</FONT>"))
                        }
                    }
                }
        }
        if (input$table == "mvNtable")
        {
                sFilename <- paste(sMarxanDir,"/output/output_mv",sep="")
                iPadding <- 5 - nchar(as.character(input$m))
                if (iPadding > 0)
                {
                        for (i in 1:iPadding)
                        {
                                sFilename <- paste(sFilename,"0",sep="")
                        }
                }
                sFilename <- paste(sFilename,input$m,".csv",sep="")
                thetable <- read.csv(sFilename,stringsAsFactors=FALSE)
                # sort the table the way spec.dat is ordered
                tableorder <- c(25,24,23,22,21,20,19,18,17,16,15,14,13,12,11,10,9,8,7,6,5,4,3,2,1)
                thetable <- thetable[tableorder,]
                # select just fields we want
                colnames(thetable)[4] <- "AmountHeld"
                colnames(thetable)[9] <- "TargetMet"
                thetable <- sqldf(paste0("SELECT Target, AmountHeld, TargetMet from thetable"))
                # load feature name from spec.dat
                specdat <- read.csv(paste0(sMarxanDir,"/input/spec.dat"),stringsAsFactors=FALSE)
                name <- sqldf(paste0("SELECT name from specdat"))
                # join the name
                thetable <- cbind(name,thetable)
                # select the fields we want from total areas file
                tafile <- read.csv(paste0(sMarxanDir,"/core1/MarOptTotalAreas.csv"))
                tafile <- tafile[tableorder,]
                tafile <- sqldf(paste0("SELECT totalarea,reservedarea from tafile"))
                # compute target gap
                targetgap <- rep(0,each=nrow(thetable))
                for (i in 1:nrow(thetable))
                {
                    if (thetable$Target[i] > 0)
                    {
                        if (thetable$AmountHeld[i] < thetable$Target[i])
                        {
                            targetgap[i] <- thetable$Target[i] - thetable$AmountHeld[i]
                        }
                    }
                }
                # join and tidy the table
                thetable <- cbind(thetable,tafile,targetgap)
                thetable <- sqldf(paste0("SELECT name, totalarea, reservedarea, Target, AmountHeld, TargetMet, targetgap from thetable"))
                colnames(thetable)[2] <- "Total"
                colnames(thetable)[3] <- "Reserved"
                colnames(thetable)[7] <- "TargetGap"
                # colour code features that have not met targets
                for (i in 1:nrow(thetable))
                {
                    if (thetable[i,6] == "no")
                    {
                        for (j in 6:7)
                        {
                            thetable[i,j] <- HTML(paste0("<FONT COLOR='blue'>",thetable[i,j],"</FONT>"))
                        }
                    } else {
                        for (j in 1:7)
                        {
                            thetable[i,j] <- HTML(paste0("<FONT COLOR='black'>",thetable[i,j],"</FONT>"))
                        }
                    }
                }
        }

        return(thetable)
    })

    output2ds <- reactive({
        input$refreshinput

        plot(sol.mds$points, xlab='', ylab='', main='NMDS of solutions', col=nmdscolours)
        text(sol.mds$points,labels=plotlabels,pos=4, col=nmdscolours)
    })

    outputdendogram <- reactive({
        input$refreshinput

        d <- dendrapply(as.dendrogram(h), labelCol)
        return(plot(d, xlab="Solutions", ylab="Disimilarity", main="Bray-Curtis dissimilarity of solutions"))
    })

    output$marxanmap <- renderPlot({
        print(outputmap())
    }, height=450,width=450)

    output$marxantable <- renderTable({
        dat <- data.frame(outputtable())
        dat
    }, sanitize.text.function = function(x) x)

    output$textfeedback = renderText({
        runmarxan()
        sprintf("Finished")
    })

    output$plot2ds <- renderPlot({ 
        print(output2ds())
    })
    
    output$plotdendogram <- renderPlot({ 
        print(outputdendogram())
    })
})
