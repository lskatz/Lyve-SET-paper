#!/bin/sh

## Usage: run this script directly passing it your target directory
## as the first argument, and it will iterate over the files in the 
## target folder and pass them as a parameter to prokkasingle.sh...

for file in $1/*
do qsub -M $2 prokkasingle.sh "$file"
done
