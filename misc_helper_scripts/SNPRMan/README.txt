Description: SNPRMan (pronounced 'Sniper-Man'): Single Nucelotide Polymorphism matrix comparisons in R via Mantel is an R script for running a Mantel test to quantify the correlation between two SNP distance matrices. This is an alternative and/or complementary method to comparing tree topologies.

Author: A. Jo Williams-Newkirk at igy7@cdc.gov
README last updated 08 Dec 2015

If you are an R wizard, please continue on to the SNPRMan documentation below. Perform whatever setup actions you like to get the listed dependencies. Be aware that to protect the muggles and throw more useful errors, SNPRMan.R requires a .Renviron file in its working directory. As an example, see example.Renviron.

If you are an R muggle (or a lazy wizard), please run the command 'Rscript setup.R /path/to/my/local/R/library' locally first. You can set up your local library in any directory with write access; make sure it exists before running setup.R. This script with download and install all R package dependencies into the directory and create a .Renviron file in your current working directory. You only have to run this script once, but SNPRMan.R will look for the .Renviron file in the working directory every time that it runs.

SNPRMan.R
Dependencies: the R packages docopt and vegan

Usage: Rscript SNPRMan.R -f <first_matrix> -s <second_matrix> [-c <corr> -p <perms> -o <outfile>]
       Rscript SNPRMan.R -v
       Rscript SNPRMan.R -h

Options:
-f <first_matrix>, --first_matrix=<first_matrix>        Tab-delimited file containing first matrix; only square matrices with both row and column labels accepted; required
-s <second_matrix>, --second_matrix=<second_matrix>     Tab-delimited file containing second matrix; only square matrices with both row and column labels accepted; required
-c <corr>, --corr=<corr>                                Correlation method either 'pearson' (parametric), 'spearman' (non-parametric), or 'kendall' (non-parametric); optional [default: spearman]
-p <perms>, --perms=<perms>                             Number of permutations; optional [default: 1000]
-o <outfile>, --outfile=<outfile>                       Output file name; optional [default: SNPRMan.output]
-h --help                                               Show this help text
-v --version                                            Show version

Files unweightedunifrac.txt and weightedunifrac.txt are included strictly as example datasets for testing. These are large matrices, so it is recommended that you set p = 10 or some other low number to avoid extended wait times for simple testing. The file cluster.submit.sh will run SNPRMan.R on the cluster using these example matrices. 
