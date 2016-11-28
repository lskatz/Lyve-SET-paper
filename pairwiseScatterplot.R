#!/usr/bin/env Rscript

library(docopt,quietly=TRUE)
commandArgs.all <- commandArgs(trailingOnly=FALSE)
commandArgs.trailingOnly <- commandArgs(trailingOnly=TRUE)

# Path to this script
thisScriptArg <- grep("--file=", commandArgs.all)
thisScript <- "SCRIPT"
if(length(thisScriptArg) > 0){
  thisScript <- substr(commandArgs.all[thisScriptArg],nchar("--file="),nchar(commandArgs.all[thisScriptArg]))
}
doc <- "Description: Compare the pairwise distances between different pipelines.

Usage: pairwiseScatterplot.R [options] <TABLE>

Where the table has at least 4 columns: isolate1 isolate2 Lyve-SETdistances otherdistances.
  Additional columns can include more 'otherdistances'

  Options:
    -h --help       Show this screen
    --max-snps=<i>  Maximum number of Lyve-SET SNPs to look at [Default: 9999]
    --image=<s>     Where to print the image [Default: pairwise.bmp]
    --limit-by=<s>  File with isolate names to include. All other isolates will be excluded.

"

opts <- docopt(doc)

library(ggplot2,quietly=TRUE)
library(reshape2,quietly=TRUE)

# Colors thanks to http://www.cookbook-r.com/Graphs/Colors_(ggplot2)
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# tall/skinny format from pairwise
t <- read.table(opts$TABLE, stringsAsFactors=FALSE, header=TRUE)
u <- reshape2::melt(t[3:length(t)], id.vars=c("Lyve.SET"), variable.name="Pipeline", value.name="Distance")
# Filter to < 50 SNPs
v <- subset(u, Lyve.SET <= as.integer(opts$"max-snps"))

# Find R2 and slope for each pipeline vs Lyve-SET
pipeline <- sort(as.character(unique(v$Pipeline)))
pipeline <- pipeline[! grepl("Lyve.SET",pipeline) ] # don't compare Lyve-SET vs Lyve-SET
cat(c("pipeline","slope","R2","\n"),sep="\t"); 
legendLabel <- c()
legendLabelStr <- ""
for (key in pipeline){
  w <- subset(v,Pipeline==key); 
  lm.tmp <- summary(lm(w$"Distance" ~ w$"Lyve.SET")); 
  slope <- lm.tmp$coefficients[2,1]; 
  r2 <- lm.tmp$r.squared; 
  intercept <- lm.tmp$coefficients[1,1]
  cat(c(key,slope,r2,"\n"),sep="\t"); 

  # Save some values. Pipeline name is still saved in `pipeline`
  legendLabel <- append(legendLabel,
    paste(key,
      " y=",signif(slope,2),"x+",signif(intercept,2), " (R2: ",signif(r2,2),")",
      sep=""
    )
  )
  #cat(c(key," y=",signif(slope,2),"x+",signif(intercept,2), " (R2: ",signif(r2,2),")","\n"),sep="")
}

legendLabelStr <- paste(legendLabel,sep="\n")

# Configure the graph
my_aes <- aes(x=Lyve.SET, y=Distance, color=factor(Pipeline))
# Make the graph
ggplot(v, my_aes) + 
  theme_minimal() + 
  #theme(legend.position="right",legend.text=element_text(size=16)) + 
  geom_point() + 
  #scale_colour_manual(values=cbPalette, label=pipeline, name="Pipeline") + 
  scale_colour_manual(values=cbPalette, label=legendLabelStr, name="Pipeline") + 
  geom_smooth(data=v, method=lm, se=FALSE, my_aes)

# Save the graph
ggsave(filename=opts$image, plot=last_plot())

