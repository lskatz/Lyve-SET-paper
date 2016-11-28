#!/bin/bash -l

# Reruns Lyve-SET from pooled vcf stage

if [ "$2" == "" ]; then
  echo "Usage: $0 Lyve-SET/msa [Lyve-SET2/msa...]"
  echo "  This script fixes any MSA that does not yet have a RAxML output in a Lyve-SET directory."
  exit 1;
fi

TMP=$(mktemp --tmpdir='.' --directory qsubSetProcessPooledVcf.XXXXXXXX)
echo "tmp dir is $TMP "

export PATH=/scicomp/groups/OID/NCEZID/DFWED/EDLB/share/bin/lyve-SET-v1.1.4e/scripts:$PATH
echo -n "Lyve-SET is being launched from ";
\which launch_set.pl

CTRL_FILE="$TMP/array.txt"

# Only keep the jobs that need the raxml output
for i in $@; do
  if [ -e "$i/out.RAxML_bipartitions" ]; then
    continue;
  fi;
  echo $i
done > $CTRL_FILE

mkdir -p $TMP/log
qsub -q all.q -N setProcessPooledVcf -o $TMP/log -j y -pe smp 1 -V -cwd -t 1-$(cat $CTRL_FILE | wc -l) \
  -v "CTRL_FILE=$CTRL_FILE" <<- "END_OF_SCRIPT"
  #!/bin/bash

  base_dir=$(sed -n ${SGE_TASK_ID}p $CTRL_FILE)
  echo "Working on $base_dir on host $(hostname)"
  mkdir -p /scratch/$USER
  scratch_out=$(mktemp --tmpdir="/scratch/$USER" --directory Lyve-SET-processPooledVcf.XXXXXX)

  cp $base_dir/out.pooled.vcf.gz $scratch_out
  tabix -f out.pooled.vcf.gz
  
  set_processPooledVcf.pl --numcpus 1 --prefix $scratch_out/out $scratch_out/out.pooled.vcf.gz
  if [ $? -gt 0 ]; then exit 1; fi;

  mv -v $scratch_out/out* $base_dir/

  rmdir $scratch_out;
END_OF_SCRIPT

