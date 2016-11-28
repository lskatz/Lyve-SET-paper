#!/bin/sh

## Embedded Grid Engine commands ##

## Designate correct shell
#$ -S /bin/sh

# The -N option sets the name of the job. This will show up in 'qstat'
#$ -N Prokka

# This sets the default directory that the script will use as its home dir
# the CWD standards for "current working diretory"
#$ -cwd

# Specifies the q
#$ -q all.q

# Send an email when the job completes
#$ -m e
##################################################################
## Now we do work ...
## From Jo's old cluster account - necessary to run.
source /etc/profile.d/modules.sh 

module load prokka/1.8

prokka "$1" --outdir "$(basename $1)" --prefix "$(basename $1)"

module unload prokka/1.8
