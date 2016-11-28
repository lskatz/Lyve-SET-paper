## Embedded Grid Engine commands ##

# Designate correct shell
#$ -S /bin/sh

# The -N option sets the name of the job.
#$ -N SNPRMan

# This sets the default directory.
#$ -cwd

# Specifies the q
#$ -q all.q

##################################################################

# Usage - For details see SNPRMan.R or run Rscript SNPRMan.R -h.
# This shell script runs the example dataset.

# For old cluster accounts.
source /etc/profile.d/modules.sh 

Rscript SNPRMan.R -f weightedunifrac.txt -s unweightedunifrac.txt -p 10
