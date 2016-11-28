#!/bin/bash -l

# Runs any TreeToReads projects in a cluster-friendly method
# Author: Lee Katz
# Usage: bash launch_TTR.sh project1 project2 [... projectX]
#   where each project has its own TTR.cfg file and associated TTR files.

if [ "$1" == "" ]; then
  echo "Usage: launch_TTR.sh project1 [project2...]"
  echo "  Where each project has its own TTR.cfg file and associated TTR files"
  exit 1;
fi

export PATH=~/bin/TreeToReads:~/bin/ART:$PATH

TMP=$(mktemp --tmpdir='.' --directory qsubTTR.XXXXXXXX)
echo "tmp dir is $TMP "

CTRL_FILE="$TMP/array.txt"
echo "$@" | tr ' ' '\n' > $CTRL_FILE

mkdir -p $TMP/log
qsub -q all.q -N TTR -o $TMP/log -j y -V -cwd -t 1-$(cat $CTRL_FILE | wc -l) \
  -v "CTRL_FILE=$CTRL_FILE" <<- "END_OF_SCRIPT"
  source /etc/profile.d/modules.sh
  module unload perl/5.16.1-MT
  export PERL5LIB=""

  export base_dir=$(sed -n ${SGE_TASK_ID}p $CTRL_FILE)
  echo "Working on $base_dir on $(hostname)"
  mkdir -p /scratch/$USER
  scratch_out=$(mktemp --tmpdir="/scratch/$USER" --directory TTR.XXXXXX)
  mkdir -p $(dirname $scratch_out)
  cd $base_dir
  sed "s|output_dir.*|output_dir = $scratch_out|" TTR.cfg > TTR.modified.cfg
  treetoreads.py TTR.modified.cfg
  if [ $? -gt 0 ]; then exit 1; fi
  cd -

  # Grab results
  mv -v $scratch_out $base_dir/TTR
  mv -v $base_dir/TTR.modified.cfg $base_dir/TTR/TTR.cfg
  
  # Shuffle the reads we so dearly wanted into 
  # the proj/reads directory.
  mkdir -p $base_dir/reads
  ls -d $base_dir/TTR/fastq/* | xargs -P $NSLOTS -n 1 sh -c 'b=$(basename $0); echo "Shuffling $b"; run_assembly_shuffleReads.pl $0/*_1.fq.gz $0/*_2.fq.gz | gzip -c > $base_dir/reads/$b.fastq.gz; echo "Finished $b"'
END_OF_SCRIPT

