#!/bin/bash

# compare all pipelines in a given directory

export dir=$1
export out="$dir/compareTrees"

if [ "$dir" == "" ]; then
  echo "This script compares the trees from a set of phylogenetic pipelines"
  echo "Usage: $0 in"
  echo "Where in is a directory, and the output folder will be in/compareTrees"
fi

export REPS=10000
NUMCPUS=12

# Include SNPRMan in the path
export PATH=$PATH:/scicomp/groups/OID/NCEZID/DFWED/EDLB/share/projects/wgsStandards/discoveringThresholds/scripts/misc_helper_scripts/SNPRMan
export PATH=$PATH:/scicomp/groups/OID/NCEZID/DFWED/EDLB/share/projects/wgsStandards/discoveringThresholds/scripts/lskScripts
source /etc/profile.d/modules.sh

if [ ! "$(which R 2> /dev/null)" ]; then
  module load R/3.2.3
fi
module load Lyve-SET/1.1.4f
module load phylip

####################################################
## Tree-level differences
####################################################
if [ -e "$out" ]; then
  echo "$out already exists! Ctrl-C or wait 2 seconds to continue."
  sleep 2
fi

mkdir -p $out

# Lyve-SET
if [ -e "$dir/Lyve-SET/msa/out.RAxML_bipartitions" ]; then
  cp -v "$dir/Lyve-SET/msa/out.RAxML_bipartitions" $out/Lyve-SET.dnd
fi

# SNP-Pipeline
if [ -e "$dir/snp-pipeline/snpma.fasta" ]; then
  cat "$dir/snp-pipeline/snpma.fasta" |\
    perl -lane '
      if(!/^>/){
        s/\-/N/g;
      }
      print;
    ' > "$dir/snp-pipeline/snpma.nogaps.fasta"
  removeUninformativeSites.pl --ambiguities-allowed "$dir/snp-pipeline/snpma.nogaps.fasta" > $dir/snp-pipeline/informative.fasta
  if [ $? -gt 0 ]; then exit 1; fi;

  if [ ! -e $dir/snp-pipeline/RAxML_bipartitions.snp-pipeline ]; then
    pushd $dir/snp-pipeline
    rm -vf RAxML*.snp-pipeline informative.fasta.phy.reduced informative.fasta.phy
    launch_raxml.sh -n $NUMCPUS informative.fasta snp-pipeline
    if [ $? -gt 0 ]; then exit 1; fi;
    popd
  fi
  
  cp -v $dir/snp-pipeline/RAxML_bipartitions.snp-pipeline $out/snp-pipeline.dnd
fi

# kSNP3
if [ -e "$dir/kSNP3/tree.majority0.75.tre" ]; then
  #KSNPTREE=$dir/kSNP3/tree.majority0.75.tre
  KSNPTREE=$dir/kSNP3/tree.ML.tre
  echo "Removing Reference from the kSNP tree at $KSNPTREE"
  pruneSafely.pl --tree $KSNPTREE reference > $out/kSNP3.dnd
  if [ $? -gt 0 ]; then
    echo "I could not prune safely (maybe 'reference' doesn't exist in $KSNPTREE); trying to copy instead";
    cp -v $KSNPTREE $out/kSNP3.dnd
  fi
  if [ $? -gt 0 ]; then echo "I could not prune safely or cp safely"; exit 1; fi;
fi

# RealPhy
if [ -e "$dir/RealPhy/reference/PolySeqOut_NoGenes" ]; then
  export REALPHYTREE="$dir/RealPhy/reference/PolySeqOut_NoGenes/polymorphisms_move.phy_phyml_tree.txt"
  echo "Removing Reference from the RealPhy tree at $REALPHYTREE";
  pruneSafely.pl --tree $REALPHYTREE reference > $out/RealPhy.dnd
  if [ $? -gt 0 ]; then echo "ERROR with RealPhy or BioPerl"; exit 1; fi;
fi

# wgMLST
if [ -e "$dir/wgMLST" ]; then
  WGMLSTTREE=$(ls $dir/wgMLST/*.dnd|head -n 1);
  cp -v $WGMLSTTREE $out/wgMLST.dnd
fi

# SNVPhyl
if [ -e "$dir/SNVPhyl" ]; then
  SNVPHYLTREE=$(ls $dir/SNVPhyl/*.newick.nhx | head -n 1);
  pruneSafely.pl --tree $SNVPHYLTREE reference > $out/SNVPhyl.dnd
  if [ $? -gt 0 ]; then
    echo "I could not prune safely (maybe 'reference' doesn't exist in $SNVPHYLTREE); trying to copy instead";
    cp -v $SNVPHYLTREE $out/SNVPhyl.dnd
    if [ $? -gt 0 ]; then echo "Could not copy $SNVPHYLTREE"; exit 1; fi;
  fi
  # if the confidence nodes are below 1, then normalize it to percentages
  confidence=$(treeInfo.pl $out/SNVPhyl.dnd --confidence|cut -f 4|tail -n 1)
  if [ "$confidence" == "" ]; then
    echo "ERROR with treeInfo.pl";
    exit 1;
  fi;
  if (( $(echo "$confidence < 1" | bc -l) )); then
    # multiply confidence by 100
    cat $out/SNVPhyl.dnd | \
    perl -MBio::TreeIO -e '
      $out=Bio::TreeIO->new(-format=>"newick"); 
      $in=Bio::TreeIO->new(-fh=>\*STDIN,-format=>"newick"); 
      while($tree=$in->next_tree){
        for my $node($tree->get_nodes){ 
          next if($node->is_Leaf); 
          next if(!$node->id);
          $node->id($node->id*100); 
        } 
      $out->write_tree($tree); 
      print "\n"; 
    }' > $out/SNVPhyl.dnd.tmp && \
    mv -v $out/SNVPhyl.dnd.tmp $out/SNVPhyl.dnd
    if [ $? -gt 0 ]; then
      echo "ERROR with converting SNVPhyl tree bootstrap values";
      exit 1;
    fi
  fi
fi

# Do some work in the output directory
# but use a subshell to ensure we don't mess with
# CWD.
(
  echo "Comparing trees with the Kendall-Colijn metric";
  ls $out/*.dnd
  cd $out
  echo 0 0.5 1 | xargs -P $NUMCPUS -n 1 sh -c '
    mkdir -pv lambda$0;
    cd lambda$0;
    /scicomp/home/gzu2/GWA/EDLB/share/projects/wgsStandards/discoveringThresholds/scripts/lskScripts/Kendall.R --rep=$REPS --rootnode=0 --lambda=$0 --plot ../Lyve-SET.dnd ../*.dnd > Kendall.tsv;
    if [ $? -gt 0 ]; then echo "ERROR with lambda=$0"; exit 1; fi;
    ls ../Lyve-SET.dnd ../*.dnd | cat -n > trees.txt
    for i in *.bmp; do
      b=$(basename $i .bmp)
      convert -verbose $i $b.png
      rm -v $i
    done
  '
  if [ $? -gt 0 ]; then exit 1; fi;
  head -n 1 lambda0/Kendall.tsv > Kendall.tsv
  sort -k1,2r -k3,3n lambda*/Kendall.tsv | grep -v flattened | grep -v Kendall | uniq >> Kendall.tsv
)
if [ $? -gt 0 ]; then exit 1; fi; # since we did a subshell

# organize just a bit
mkdir -v $out/flattenedTrees
cp -v $out/lambda0/*.dnd $out/flattenedTrees/
rm -v $out/lambda*/*.dnd

echo "Running treedist with both Kuhner-Felstein and Robinson-Foulds"
for method in kf rf; do
  statsOutFile=$out/biophylo.$method.tsv;
  rm -vf $statsOutFile
  /scicomp/home/gzu2/GWA/EDLB/share/projects/wgsStandards/discoveringThresholds/scripts/lskScripts/treedist_wrapper.pl --method $method --numtrees $REPS --numcpus $NUMCPUS $out/Lyve-SET.dnd $out/*.dnd > $statsOutFile
  if [ $? -gt 0 ]; then 
    echo "ERROR with treedist_wrapper.pl with method $method and $out/*.dnd"; 
    rm -vf $statsOutFile
    continue;
  fi;
  sed -i "s|$out/||g" $statsOutFile;
  if [ $? -gt 0 ]; then 
    echo "ERROR with sed on $statsOutFile";
    continue;
  fi;
done

