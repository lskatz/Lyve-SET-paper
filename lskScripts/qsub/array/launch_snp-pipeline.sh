#!/bin/bash -l

# Runs any set of reads through Snp-Pipeline in a cluster-friendly way.
# Each reads directory will be a distinct project.
# Author: Lee Katz
# Usage: bash $0 reference.fasta readsdir1 readsdir2 [... readsdirX]

source /etc/profile.d/modules.sh
#module load snp-pipeline/0.5.2
source ~/.local/virtualenv/python2.7/bin/activate  # update python and include local snp-pipeline
module load java/latest
module load perl/5.16.1-MT
export PATH=$(sed 's/Lyve-SET//i' <<< $PATH);  # just jack up any Lyve-SET path to avoid mergeVcf namespace conflicts
export PATH=$PATH:/scicomp/groups/OID/NCEZID/DFWED/EDLB/share/bin/CG-Pipeline/scripts

if [ "$2" == "" ]; then
  echo "Usage: $0 ref/reference.fasta dir [dir2 ... ]"
  exit 1;
fi

TMP=$(mktemp --tmpdir='.' --directory qsubSnp-Pipeline.XXXXXXXX)
echo "tmp dir is $TMP "

REF=$1; shift; # get the reference genome and remove it from ARGV

CTRL_FILE="$TMP/array.txt"
echo "$@" | tr ' ' '\n' | grep . > $CTRL_FILE

mkdir -p $TMP/log
# Has to have an exclusive node because it is a greedy script
#qsub -q all.q -N snp-pipeline -o $TMP/log -j y -pe smp 12,16 -hard -l exclusive=true -V -cwd -t 1-$(cat $CTRL_FILE | wc -l) \
qsub -q all.q -N snp-pipeline -o $TMP/log -j y -pe smp 12,16 -V -cwd -t 1-$(cat $CTRL_FILE | wc -l) \
  -v "CTRL_FILE=$CTRL_FILE" -v "REF=$REF" <<- "END_OF_SCRIPT"
  #!/bin/bash

  module load snp-pipeline/0.5.2
  module load java/latest
  module load perl/5.16.1-MT

  base_dir=$(sed -n ${SGE_TASK_ID}p $CTRL_FILE)
  echo "Working on $base_dir on $(hostname)"
  mkdir -p /scratch/$USER
  scratch_out=$(mktemp --tmpdir="/scratch/$USER" --directory snp-pipeline.XXXXXX)
  export scratch_out
  mkdir -p $scratch_out/samples

  # Make the config file and put it into scratch_out (snppipeline.conf)
  copy_snppipeline_data.py configurationFile $scratch_out
  if [ $? -gt 0 ]; then echo "ERROR copying configuration file"; exit 1; fi;

  # Find what reads we're using. Assume all reads have been shuffled.
  READS=$(ls $base_dir/reads/*.fastq.gz);

  # Deshuffle reads into the correct directories
  echo $READS | xargs -P $NSLOTS -n 1 sh -c '
    sample=$(basename $0 .fastq.gz);
    sampleDir=$scratch_out/samples/$sample;
    mkdir -p $sampleDir;
    echo "Deshuffling into $sampleDir";
    run_assembly_shuffleReads.pl $0 -d -gz 1>$sampleDir/1.fastq.gz 2>$sampleDir/2.fastq.gz
    if [ $? -gt 0 ]; then echo "ERROR deshuffling $0"; exit 1; fi;
  ';

  # Run snp-pipeline on the scratch drive without qsub inception.
  run_snp_pipeline.sh -c $scratch_out/snppipeline.conf -s $scratch_out/samples -m copy -o $scratch_out $REF
  if [ $? -gt 0 ]; then echo "ERROR with run_snp_pipeline.sh"; exit 1; fi;
  rm -rvf $scratch_out/samples/*

  rm -rfv $base_dir/snp-pipeline
  mv -v $scratch_out $base_dir/snp-pipeline
END_OF_SCRIPT

