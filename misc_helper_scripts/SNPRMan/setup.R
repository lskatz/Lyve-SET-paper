#!/usr/bin/Rscript

# Initial set up script for SNPRMan. Run LOCALLY to install all required packages in your local library BEFORE running SNPRMan. Requires a single positional argument which is the path to your local library.
# Author: A. Jo Williams-Newkirk at igy7@cdc.gov
# Last updated: 07 December 2015
# Version 0.1

# Parse args
args <- commandArgs(trailingOnly = TRUE)
if(length(args) != 1) {
  stop("Setup requires a single positional argument, which is the path to your local R library. Please try again.")
}

# Check to see if library path exists
if(!file.exists(args[1])) {
    stop("You have entered a library path that does not exist. Please create the folder before running Setup or check for typos.")
}

# Install all dependencies to local library if not already installed.
print("This may take a while if you are missing many of the required packages.")
pks <- c("permute", "lattice", "MASS", "cluster", "mgcv", "vegan", "stringi", "magrittr", "stringr", "docopt")
for(p in 1:length(pks)) {
  if(!pks[p] %in% installed.packages(lib.loc = args[1])) {
    install.packages(pkgs = pks[p], 
                     lib = args[1],
                     repos = "http://cran.us.r-project.org",
                     verbose = FALSE,
                     quiet = TRUE)
    print(paste(pks[p], "was successfully installed."))
  } else {
      print(paste(pks[p], "was already installed. Skipping."))
  }
}
print("Package installation complete.")

# Check to see if a .Renviron file already exists in the current working directory; create if it does not.
if(file.exists(".Renviron")) {
  stop("Setup detects an existing .Renviron file in your current working directory. To avoid corrupting an existing file, Setup suggests that you manually inspect this file and add 'R_LIBS=' followed by the path to your local R library if it is not already present.")
} else {
  fileConn <- file(".Renviron")
  writeLines(paste("R_LIBS=", args[1], sep = ""), 
             fileConn)
  close(fileConn)
  print("Created .Renviron file.")
}

print("SNPRMan Setup is complete. You should only have to run Setup once as long as you keep a copy of .Renviron in your working directory. Have a nice day.")
