#!/usr/bin/env Rscript

library(docopt,quietly=TRUE)
commandArgs.all <- commandArgs(trailingOnly=FALSE)
commandArgs.trailingOnly <- commandArgs(trailingOnly=TRUE)

doc <- "Description: Compare the pairwise distances between different pipelines.

Usage: pairwiseScatterplot.R [options] <TABLE>...

Where the table has at least 4 columns: isolate1 isolate2 Lyve-SET otherdistances.
  Additional columns can include more 'otherdistances'

  Options:
    -h --help       Show this screen
    --max-snps=<i>  Maximum number of Lyve-SET SNPs to look at [Default: 9999]
    --max-dist=<i>  Maximum number of other-pipeline SNPs to look at [Default: 9999]
    --image=<s>     Where to print the image [Default: pairwise.bmp]
    --limit-by=<s>  File with isolate names to include. All other isolates will be excluded.

"

opts <- docopt(doc)

library(ggplot2,quietly=TRUE)
library(reshape2,quietly=TRUE)
library(plyr, quietly=TRUE)
library(gridExtra, quietly=TRUE)

logmsg <- function(...){
  # Init
  defaultScriptName <- "SCRIPT.R"
  thisScript <- defaultScriptName  # to be changed later, hopefully
  
  # Check if it's even being run via Rscript with arguments
  if(!exists("commandArgs.all")){
    commandArgs.all <- "--file="+defaultScriptName
  }
  
  # Get the path to this script
  thisScriptArg <- grep("--file=", commandArgs.all)
  if(length(thisScriptArg) > 0){
    # For the argument that has --file=, substring away the --file= part.
    thisScript <- substr(commandArgs.all[thisScriptArg],nchar("--file=")+1,nchar(commandArgs.all[thisScriptArg]))
    # Remove the directory name
    thisScript <- basename(thisScript)
  }

  cat(thisScript,": ",..., "\n" ,sep="", file=stderr())
}

# http://stackoverflow.com/questions/2104483/how-to-read-table-multiple-files-into-a-single-table-in-r
read.tables <- function(file.names, ...) {
  ldply(file.names, function(fn) data.frame(Filename=fn, read.csv(fn, ...)))
}

# Colors thanks to http://www.cookbook-r.com/Graphs/Colors_(ggplot2)
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

t <- read.tables(opts$TABLE, sep="\t")
t$Filename <- sapply(t$Filename,FUN=function(filename){
                     filename <- basename(dirname(as.character(filename)));
                     filename <- gsub('LinearRegression','',filename,perl=TRUE)
                     return(filename)
              })
t <- t[c(1,4:length(t))]   # remove names of genomes

u <- reshape2::melt(t, id.vars=c("Filename","Lyve.SET"), variable.name="Pipeline", value.name="Distance")

# Disregard Lyve-SET vs Lyve-SET if it exists
v <- subset(u, Pipeline != "Lyve.SET")
# Filter to < 50 SNPs
v <- subset(v, Lyve.SET <= as.integer(opts$"max-snps"))
v <- subset(v, Distance <= as.integer(opts$"max-dist"))

# Make a four-panel graph
graphVector <- c()
for (filename in as.character(unique(v$Filename))){
  # Make a legend label and calculate stope vs Lyve-SET
  legendLabel <- c()
  for (pipeline in as.character(unique(v$Pipeline))){
    w <- subset(v,Pipeline==pipeline); 
    w <- subset(w,Filename==filename)

    error <- 0;  # assume no error until told otherwise
    tryCatch({
      lm.tmp <- summary(lm(w$"Distance" ~ w$"Lyve.SET")); 
    }, error = function(e){
      errstr <- gsub('^\\s+|\\s+$','',as.character(e), perl=TRUE)  # whitespace trim
      #logmsg("Could not get a graph for ",filename,"/",pipeline,": ",errstr);
      error <<- 1
    })

    if(error > 0){
      next;
    }

    slope <- lm.tmp$coefficients[2,1]; 
    r2 <- lm.tmp$r.squared; 
    intercept <- lm.tmp$coefficients[1,1]
    cat(c(filename,pipeline,slope,r2,"\n")); 


    # Save some values. Pipeline name is still saved in `pipeline`
    legendLabel <- append(legendLabel,
      paste(pipeline,
        " y=",signif(slope,2),"x+",signif(intercept,2), " (R2: ",signif(r2,2),")",
        sep=""
      )
    )
  }
  legendLabelStr <- paste(legendLabel,sep="\n")

  my_aes <- aes(x=Lyve.SET, y=Distance, color=factor(Pipeline))
  my_graph <- ggplot(v, my_aes) + 
    theme_minimal() + 
    theme() + 
    geom_point() + 
    scale_colour_manual(values=cbPalette, label=legendLabelStr, name="Pipeline") + 
    geom_smooth(data=v, method=lm, se=FALSE, my_aes)

  graphVector <- append(graphVector,my_graph)
}

my_grid <- do.call(grid.arrange, c(graphVector[1],graphVector[2], graphVector[3], graphVector[4], ncol=2, top="Main title"))
#ggsave(filename=opts$image, plot=my_grid)
# Make a grid
#my_graph <- my_graph + facet_wrap(~ Filename, ncol=2, nrow=2)
#ggsave(filename=opts$image, plot=my_graph)



