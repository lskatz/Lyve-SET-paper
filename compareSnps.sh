#!/bin/bash

# compare all pipelines in a given directory

dir=$1
out="$dir/compareSnps";

if [ "$dir" == "" ]; then
  echo "This script compares the trees from a set of phylogenetic pipelines"
  echo "Usage: $0 in"
  echo "Where in is a directory, and the output folder will be in/compareTrees"
fi

REPS=10000
NUMCPUS=12

mkdir -p $out/positions
# Include SNPRMan in the path
export PATH=$PATH:/scicomp/groups/OID/NCEZID/DFWED/EDLB/share/projects/wgsStandards/discoveringThresholds/scripts/misc_helper_scripts/SNPRMan:/scicomp/groups/OID/NCEZID/DFWED/EDLB/share/projects/wgsStandards/discoveringThresholds/scripts
export PATH=/scicomp/groups/OID/NCEZID/DFWED/EDLB/share/projects/wgsStandards/discoveringThresholds/scripts/lskScripts:$PATH
# Lyve-SET
export PATH=/scicomp/groups/OID/NCEZID/DFWED/EDLB/share/bin/lyve-SET-v1.1.4f/scripts:$PATH
source /etc/profile.d/modules.sh
module load perl/5.16.1-MT
module load tabix
if [ ! "$(which R 2> /dev/null)" ]; then
  module load R/3.2.3
fi

############################################
## Now for SNP-level differences
############################################

TSV=""
# Lyve-SET
if [ -e "$dir/Lyve-SET/msa/out.snpmatrix.tsv" ]; then 
  if [ ! -e $out/Lyve-SET.tsv ]; then
    filterMatrix.pl "$dir/Lyve-SET/msa/out.snpmatrix.tsv" --ambiguities --Ns-as-ref > $out/positions/Lyve-SET.tsv
    if [ $? -gt 0 ]; then echo "Problem with filterMatrix.pl"; exit 1; fi;
    cat $out/positions/Lyve-SET.tsv | matrixToAlignment.pl | pairwiseDistances.pl -n $NUMCPUS | sort -k3,3n | pairwiseTo2d.pl | perl -lane 's/\t-\t/\t0\t/g; print;' > $out/Lyve-SET.tsv
    if [ $? -gt 0 ]; then echo "ERROR converting Lyve-SET positions to distances"; exit 1; fi;
  fi
  TSV="$TSV $out/Lyve-SET.tsv"
fi

if [ -e "$dir/snp-pipeline/snpma.fasta" ]; then
  if [ ! -e $out/snp-pipeline.tsv ]; then
    cut -f 1,2 "$dir/snp-pipeline/snplist.txt" > $out/positions/snp-pipeline.tsv
    if [ $? -gt 0 ]; then "ERROR getting snp-pipeline positions"; exit 1; fi;
    pairwiseDistances.pl -n $NUMCPUS $dir/snp-pipeline/snpma.fasta | sort -k3,3n  | pairwiseTo2d.pl | perl -lane 's/\t-\t/\t0\t/g; print;' > $out/snp-pipeline.tsv
    if [ $? -gt 0 ]; then exit 1; fi;
  fi
  TSV="$TSV $out/snp-pipeline.tsv"
fi

#if [ -e "$dir/kSNP3/SNPs_in_majority0.75_matrix.fasta" ]; then
if [ -e "$dir/kSNP3/SNPs_all_matrix.fasta" ]; then
  if [ ! -e $out/kSNP3.tsv ]; then
    KSNPVCF=$(ls $dir/kSNP3/VCF.*.vcf | head -n 1)
    ksnpsToVcf.pl $dir/kSNP3/SNPs_all > $dir/kSNP3/kSNP3.vcf
    if [ $? -gt 0 ]; then echo "ERROR with creating a kSNP3 VCF"; exit 1; fi;
    fixKsnpVcf.pl -ref $dir/../reference/reference.fasta < $dir/kSNP3/kSNP3.vcf | vcf-sort | bgzip -c > $out/positions/kSNP3.vcf.gz && \
    tabix $out/positions/kSNP3.vcf.gz
    if [ $? -gt 0 ]; then echo "ERROR with fixing kSNP3 VCF"; exit 1; fi;

    pooledToMatrix.sh -o $out/positions/kSNP3.tsv $out/positions/kSNP3.vcf.gz
    if [ $? -gt 0 ]; then echo "ERROR with making a kSNP3 snp matrix"; exit 1; fi;
    #pairwiseDistances.pl -n $NUMCPUS $dir/kSNP3/SNPs_in_majority0.75_matrix.fasta | grep -v 'reference' | sort -k3,3n  | pairwiseTo2d.pl | perl -lane 's/\t-\t/\t0\t/g; print;' > $out/kSNP3.tsv
    pairwiseDistances.pl -n $NUMCPUS $dir/kSNP3/SNPs_all_matrix.fasta | grep -v 'reference' | sort -k3,3n | pairwiseTo2d.pl | perl -lane 's/\t-\t/\t0\t/g; print;' > $out/kSNP3.tsv
    if [ $? -gt 0 ]; then exit 1; fi;
  fi
  TSV="$TSV $out/kSNP3.tsv"
fi

if [ -e "$dir/RealPhy/reference/PolySeqOut_NoGenes/polymorphisms_move.fas" ]; then
  if [ ! -e $out/RealPhy.tsv ]; then
    # custom perl to make a VCF file
    for i in $dir/RealPhy/reference/PolySeqOut_NoGenes/*details.txt; do
      export b=$(basename $i details.txt);
      tail -n +2 $i | perl -lane 'BEGIN{ print "##fileformat=VCFv4.0\n##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n##reference=$dir/../reference/reference.fasta\n". join("\t","#CHROM","POS","ID","REF","ALT","QUAL", "FILTER","INFO","FORMAT",$ENV{b});} print join("\t",$F[1],$F[4],".",$F[6],$F[7],".","PASS",".","GT","1")' | bgzip -c > $dir/RealPhy/reference/PolySeqOut_NoGenes/$b.vcf.gz
      tabix $dir/RealPhy/reference/PolySeqOut_NoGenes/$b.vcf.gz
    done
    mergeVcf.sh -o $out/positions/RealPhy.vcf.gz $dir/RealPhy/reference/PolySeqOut_NoGenes/*.vcf.gz
    pooledToMatrix.sh -o $out/positions/RealPhy.tsv $out/positions/RealPhy.vcf.gz
    if [ $? -gt 0 ]; then echo "ERROR with pooledToMatrix.sh"; exit 1; fi;

    pairwiseDistances.pl -n $NUMCPUS $dir/RealPhy/reference/PolySeqOut_NoGenes/polymorphisms_move.fas | grep -v "reference" | sort -k3,3n  | pairwiseTo2d.pl | perl -lane 's/\t-\t/\t0\t/g; print;' > $out/RealPhy.tsv
    if [ $? -gt 0 ]; then exit 1; fi;
  fi
  TSV="$TSV $out/RealPhy.tsv"
fi

if [ -e "$dir/SNVPhyl" ]; then
  # SNP positions
  thispositions=$(ls $dir/SNVPhyl/*pseudo-positions.tsv.tabular|head -n 1)
  if [ ! -e "$thispositions" ]; then
    echo "ERROR: cannot find SNVPhyl snp positions at $dir/SNVPhyl/*pseudo-positions.tsv.tabular"
    exit 1
  fi

  cat $thispositions | perl -MList::Util=uniq -lane 'BEGIN{ $header=<>; chomp($header); @sample=split(/\t/,$header); splice(@sample,0,4); print "##fileformat=VCFv4.0\n##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n##reference=$dir/../reference/reference.fasta\n". join("\t","#CHROM","POS","ID","REF","ALT","QUAL", "FILTER","INFO","FORMAT",@sample);} next if($F[2] ne "valid"); @GT=@F[4..@F-1]; @ALT=grep {$_ ne $F[3]} uniq(@GT); $ALT=join(",",@ALT); print join("\t",$F[0],$F[1],".",$F[3],$ALT,".", ".", ".", ".", @GT);' | bgzip -c > $out/positions/SNVPhy.vcf.gz
  tabix $out/positions/SNVPhy.vcf.gz
  if [ $? -gt 0 ]; then exit 1; fi;

  # SNP distances
  thistsv=$(ls $dir/SNVPhyl/*snp_matrix.tsv.tabular | head -n 1)
  if [ ! -e "$thistsv" ]; then
    echo "ERROR: cannot find SNVPhyl snp matrix at $dir/SNVPhyl/*snp_matrix.tsv.tabular"
    exit 1
  fi
  # Remove the column and row with "reference"
  index=$(head -n 1 $thistsv | perl -lane 'for(my $i=0;$i<@F;$i++){if($F[$i] eq "reference"){print $i+1; exit 0;}} print "not found"')
  cut -f 1-$(($index-1)),$(($index+1))- $thistsv | sed 's/strain/\./' | grep -v reference > $out/SNVPhyl.tsv
  TSV="$TSV $out/SNVPhyl.tsv"
fi

# wgMLST
if [ -e "$dir/wgMLST" ]; then
  #thistsv=$(ls $dir/wgMLST/*similaritymatrix.txt | head -n 1)
  alleles=$(ls $dir/wgMLST/*alleleMatrix.csv|head -n 1);
  if [ ! -e "$alleles" ]; then
    echo "ERROR: cannot find the alleles file to make pairwise distances at $dir/wgMLST/*alleleMatrix.csv";
    exit 1
  fi

  pairwiseAlleleDistances.pl $alleles > $dir/wgMLST/matrix.tsv && cp -v $dir/wgMLST/matrix.tsv $out/wgMLST.tsv
  if [ $? -gt 0 ]; then exit 1; fi;

  TSV="$TSV $out/wgMLST.tsv"

  #cat $thistsv | perl -MData::Dumper -Mautodie -MBio::Matrix::IO -e 'open(TRANSL,"key_translation.csv"); while(<TRANSL>){ next if(/Key/ && /ID/); s/^\s+|\s+$//g; ; @F=split /,/; $transl{$F[2]}=$F[9]; } $in=Bio::Matrix::IO->new(-fh=>\*STDIN, -format=>"phylip", -verbose=>0)->next_matrix; @samples=sort {$a cmp $b} keys(%transl); $numSamples=@samples; print join("\t",".",map { $transl{$_} }@samples)."\n"; for(my $i=0;$i<$numSamples;$i++){ print $transl{$samples[$i]}; for( my $j=0;$j<$numSamples;$j++){ print "\t".($in->get_entry($samples[$i],$samples[$j]) || $in->get_entry($samples[$j],$samples[$i]) );} print "\n"; } ' > $out/wgMLST.tsv
  #TSV="$TSV $out/wgMLST.tsv";
fi

# Run the Mantel test between all results in the output directory
(
for i in $TSV; do
  export i;
  export out;
  export REPS
  numSamples=$(($(wc -l < $i) - 1))
  for j in $TSV; do
    echo "Running SNIPRMAN on $i and $j";
    if [ "$i" == "$j" ]; then continue; fi;

    # sort the matrix to avoid the sort bug in SNPRMan
    for matrix in $i $j; do
      (echo "    $numSamples"; tail -n +2 $matrix) | perl -MData::Dumper -MBio::Matrix::IO -e '$matrix=Bio::Matrix::IO->new(-fh=>\*STDIN, -format=>"phylip")->next_matrix; @names=grep { $_ ne "." } sort {$a cmp $b} @{ $matrix->names }; print join("\t",".",@names)."\n"; for(my $i=0;$i<@names;$i++){ print $names[$i]; for(my $j=0;$j<@names;$j++){print "\t".$matrix->get_entry($names[$i],$names[$j]); } print "\n"; }'  > $matrix.tmp
    done

    outfile=$(basename $i .tsv)_$(basename $j .tsv).snprman;
    SNPRMan.R -f $i.tmp -s $j.tmp -p $REPS -o $out/$outfile
    if [ $? -gt 0 ]; then exit 1; fi;
    rm $i.tmp $j.tmp

  done

  echo "Creating a matrix image for $i"
  distanceMatrixToImage.pl --px-per-square 16 --left-margin 200 --sort-by-matrix $out/Lyve-SET.tsv --max-distance 100 --outfile $i.100max.bmp $i
  if [ $? -gt 0 ]; then exit 1; fi;
done
)
if [ $? -gt 0 ]; then exit 1; fi;

# Find the r^2 values
grep 'r^2' $out/*.snprman | perl -lane 'tr/:_/\t\t/; s/\.out|r\^2\t//g; print;' | sort -k1,1 -k2,2 > $out/snprman.tsv &&\
rm -vf $out/*.snprman

pipelineSnpsVsPipeline.pl --lyveset $dir/Lyve-SET --snppipeline $dir/snp-pipeline --RealPhy $dir/RealPhy --ksnp3 $dir/kSNP3 --snvphyl $dir/SNVPhyl --reference $dir/../reference/reference.fasta > $out/positions.tsv
if [ $? -gt 0 ]; then
  exit 1;
fi;

# Make a scatterplot of pairwise distances
# Make tall/skinny tables
listOfOutbreak=/scicomp/groups/OID/NCEZID/DFWED/EDLB/share/projects/wgsStandards/discoveringThresholds/coverage/combinedOutbreaks/outbreakIsolates.tsv
for i in $TSV; do 
  if [ ! -e "$i" ]; then continue; fi; 
  b=$(basename $i .tsv); 
  outfile=$out/$b.tallskinny.tsv; 
  outbreakOutfile="$out/$b.outbreakOnly.tallskinny.tsv"
  perl -Mautodie -lane '
    BEGIN{
      # Treat the first line of the matrix differently
      $header=<>; 
      chomp($header); 
      @header=split(/\t/,$header); 
      shift(@header); 
    } 
    $thisName=shift(@F); 
    for(my $i=0;$i<@F;$i++){ 
      my($n1,$n2)=sort {$a cmp $b} ($thisName, $header[$i]); 
      next if($n1 eq $n2); 
      $F[$i]=0 if($F[$i] eq "-"); 
      print join("\t",$n1,$n2,$F[$i]); 
  }' < $i | sort | uniq > $outfile; 

  perl -Mautodie -lane '
    BEGIN{
      # Read the outbreak isolates file
      open(IN,"'$listOfOutbreak'"); 
      @outbreak=<IN>; 
      close IN; 
      chomp(@outbreak);
      @outbreak{@outbreak}=(1) x scalar(@outbreak);
      
      $header=<>;
      @header=split /\t/, $header;
      $numHeaders=@header;
      chomp($header);
      print $header;
    }
    if(!$outbreak{$F[0]}){
      #print STDERR "$F[0] not found";
      next;
    }elsif(!$outbreak{$F[1]}){
      #print STDERR "$F[1] not found";
      next;
    }

    print;
  ' < $outfile > $outbreakOutfile
  
done;

# Combine the tall/skinny tables
(echo -e "g1\tg2\tLyve-SET"; sort -k1,2 $out/Lyve-SET.tallskinny.tsv) > $out/tmp.tsv; for i in $TSV; do b=$(basename $i .tsv); if [ $b == "Lyve-SET" ]; then continue; fi; tall=$out/$b.tallskinny.tsv; if [ ! -e $tall ]; then echo "Not found: $tall"; continue; fi; paste $out/tmp.tsv <(echo $b; sort -k1,2 $tall | cut -f 3) > $out/tmp.tsv.tmp && mv $out/tmp.tsv.tmp $out/tmp.tsv; done;
#TODO Add in outbreak or nonoutbreak
mv $out/tmp.tsv $out/pairwiseCorrelation.tsv

(
  # Make the outbreaks-only file
  perl -Mautodie -lane '
    BEGIN{
      open(IN,"'$listOfOutbreak'");
      @outbreak=<IN>; 
      close IN; 
      chomp(@outbreak); 
      @outbreak{@outbreak}=(1) x scalar(@outbreak); 
      chomp($header=<>); 
      print $header; 
      print STDERR scalar(@outbreak) ." isolates in the outbreak table";

      $header=<>;
      @header=split /\t/, $header;
      $numHeaders=@header;
      chomp($header);
      print $header;
    } 
    next if(!$outbreak{$F[0]} || !$outbreak{$F[1]} || @F<$numHeaders); 
    
    # Stay under the Jackson et al recommendation?
    $tightDistance=1;
    for(2..$numHeaders){
      $tightDistance=0 if($F[$_] > 100);
    }
    next if(!$tightDistance);

    $num++; 
    print; 
    END{
      print STDERR "$num pairwise distances printed";
    }
  ' < $out/pairwiseCorrelation.tsv > $out/pairwise.outbreakOnly.tsv

  cd $out
  Rscript ../../../../scripts/pairwiseScatterplot.R --image pairwiseCorrelation.png  pairwiseCorrelation.tsv && \
  Rscript ../../../../scripts/pairwiseScatterplot.R --max-snps 50 --image pairwise.outbreakOnly.png pairwise.outbreakOnly.tsv

  if [ $? -gt 0 ]; then
    echo "ERROR with pairwiseScatterplot.R"
    exit 1;
  fi
)



