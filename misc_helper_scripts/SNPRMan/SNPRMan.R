#!/usr/bin/Rscript

# Last updated: 08 December 2015
# Version
ver = "0.1"

### Functions

# Import matrix file
impDist <- function(infile) {
  # read in file
  tmp <- read.delim(file = infile, 
                    header = TRUE, 
                    row.names = 1, 
                    sep = "\t", 
                    stringsAsFactors = FALSE)
  # check for and remove empty columns from badly formatted files
  if(all(is.na(tmp[,length(colnames(tmp))])) == TRUE) {
    tmp[,length(colnames(tmp))] <- NULL
  }
  # change any '-' to proper 0's
  tmp[tmp == '-'] = '0'
  # change data type to numeric
  tmp <- apply(tmp, 2, as.numeric)
  # fix the row names
  rownames(tmp) <- colnames(tmp)
  tmp <- as.dist(tmp)
  return(tmp)
}

### End functions

# Load required libraries.
suppressMessages(library(vegan))  # Mantel function
suppressMessages(library(docopt)) # arg parsing

# Configure arg parsing
doc <- "Description: SNPRMan (pronounced 'Sniper-Man'): Single Nucelotide Polymorphism matrix comparisons in R via Mantel is an R script for running a Mantel test to quantify the correlation between two SNP distance matrices. This is an alternative and/or complementary method to comparing tree topologies. See the included README for addional set up instructions.

Author: A. Jo Williams-Newkirk at igy7@cdc.gov.

Dependencies:
  R packages: docopt, vegan

Usage: SNPRMan.R -f <first_matrix> -s <second_matrix> [-c <corr> -p <perms> -o <outfile>]
       SNPRMan.R -v
       SNPRMan.R -h

Options:
-f <first_matrix>, --first_matrix=<first_matrix>        Tab-delimited file containing first matrix; only square matrices with both row and column labels accepted; required
-s <second_matrix>, --second_matrix=<second_matrix>     Tab-delimited file containing second matrix; only square matrices with both row and column labels accepted; required
-c <corr>, --corr=<corr>                                Correlation method either 'pearson' (parametric), 'spearman' (non-parametric), or 'kendall' (non-parametric); optional [default: spearman]
-p <perms>, --perms=<perms>                             Number of permutations; optional [default: 1000]
-o <outfile>, --outfile=<outfile>                       Output file name; optional [default: SNPRMan.output]
-h --help                                               Show this help text
-v --version                                            Show version"

# Parse args
opt <- docopt(doc = doc, version = ver)

# Import matrix files.
m1 <- impDist(opt$first_matrix)
m2 <- impDist(opt$second_matrix)

# Run Mantel test
test <- mantel(xdis = m1, 
        ydis = m2, 
        method = opt$corr, 
        permutations = as.numeric(opt$perms))

# Write log file
fileConn <- file(opt$outfile)
writeLines(c(paste("SNPRMan.R. version", ver),
             R.version.string, 
             paste("Package Vegan version:", packageVersion("vegan")),
             paste("Time:", Sys.time()),
             paste("Date:", Sys.Date()),
             paste("First matrix file:", opt$first_matrix),
             paste("Second matrix file:", opt$second_matrix),
             paste("Correlation method requested:", opt$corr),
             paste("Correlation method used:", test$method),
             paste("Number of permutations requested:", opt$perms),
             paste("Number of permutations performed:", test$permutations),
             paste("Mantel statistic r:", test$statistic),
             paste("r^2:", as.numeric(test$statistic)^2),
             paste("P-value (by permutation):", test$signif)
             ),
           fileConn)
close(fileConn)
