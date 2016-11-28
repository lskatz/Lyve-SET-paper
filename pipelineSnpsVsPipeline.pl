#!/usr/bin/env perl
# Compares sets of SNPs
#

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use File::Basename qw/basename/;

local $0=basename $0;
sub logmsg{print STDERR "$0: @_\n"}
exit main();

sub main{
  my $settings={};
  GetOptions($settings,qw(help numcpus=i ref|reference=s lyveset=s@ ttr=s@ snppipeline=s@ realphy=s@ ksnp3=s@ snvphyl=s@)) or die $!;
  $$settings{numcpus}||=1;
  die usage() if($$settings{help});

  my %project;
  for my $type(qw(lyveset ttr snppipeline realphy ksnp3 snvphyl)){
    # Copy the array over instead of by reference
    $$settings{$type}//=[];
    $project{$type}=[ @{ $$settings{$type} } ];
  }
  my @project=flatten(values(%project));
  my $numProjects=@project;

  if($numProjects < 2){
    die "ERROR: need at least two projects for comparison\n".usage();
  }

  # Determine positions of SNP sites and maybe their
  # allele calls.
  my %snpSites; # project=>{site1=>A,site2,...}
  while(my($projectType,$projectList)=each(%project)){
    for my $project(@$projectList){
      $snpSites{$project}=getSnpSites($projectType,$project,$settings);
    }
  }

  print join("\t",qw(ref query TP TN FP FN Sn Sp))."\n";
  for(my $i=0;$i<$numProjects;$i++){
    next if(!keys(%{ $snpSites{$project[$i]} }));
    for(my $j=0; $j<$numProjects; $j++){
      next if($i==$j); # no need for self vs self

      next if(!keys(%{ $snpSites{$project[$j]} }));

      my($TP,$TN,$FP,$FN,$Sn,$Sp) = compareSnps($snpSites{$project[$i]},$snpSites{$project[$j]},$settings);
      my $ref=basename($project[$i]);
      my $query=basename($project[$j]);
      $_=sprintf("%0.2f",$_) for($Sn, $Sp);
      print join("\t",$ref,$query,$TP,$TN,$FP,$FN,$Sn,$Sp)."\n";
    }
  }


  return 0;
}

sub getSnpSites{
  my($projectType,$project,$settings)=@_;
  my $snpSites={};
  if($projectType eq 'lyveset'){
    $snpSites=lyveSetSnps($project,$settings);
  } elsif($projectType eq 'ttr'){
    $snpSites=ttrSnps($project,$settings);
  } elsif($projectType eq 'snppipeline'){
    $snpSites=snppipelineSnps($project,$settings);
  } elsif($projectType eq 'ksnp3'){
    $snpSites=ksnp3Snps($project,$settings);
  } elsif($projectType eq 'snvphyl'){
    $snpSites=snvphylSnps($project,$settings);
  } elsif($projectType eq 'realphy'){
    $snpSites=realphySnps($project,$settings);
  } else {
    die "ERROR: I do not understand project type $projectType";
  }

  return $snpSites;
}

sub compareSnps{
  my ($realSnps,$querySnps,$settings)=@_;

  # Initialize counts to zero
  my($TP,$TN,$FP,$FN)=(0) x 4;

  my $numSites=keys(%$realSnps);
  
  # How many of the real SNPs were found?
  while(my($trueChrom,$truePosHash)=each(%$realSnps)){
    while(my($truePos,$trueRef)=each(%$truePosHash)){
      # If the site is marked as a nonsnp and if the query
      # does not report it, then it is a true negative.
      if(!$trueRef && !$$querySnps{$trueChrom}{$truePos}){
        $TN++;
      } elsif($$querySnps{$trueChrom}{$truePos}){
        $TP++;  # True positive: correct SNP was found
      } else {
        $FN++;  # False negative: a real SNP was not found
      }
    }
  }

  # How many SNPs were found that were not real?
  while(my($chrom,$queryRefHash)=each(%$querySnps)){
    while(my($queryPos,$queryRef)=each(%$queryRefHash)){
      # For a position found in the query,
      # if it is a SNP and it is not found in
      # the real SNPs, then it is a FP
      if($queryRef && !$$realSnps{$chrom}{$queryPos}){
        $FP++;
      }
    }
  }

  # Sensitivity and specificity
  die "ERROR: there are no true positives and no false negatives. This can happen if your results did not use the same reference genome" if($TP + $FN < 1);
  my $Sn=$TP/($TP+$FN);
  my $Sp=0;
  if($TN+$FP > 0){
    $Sp=$TN/($TN+$FP);
  }

  return ($TP,$TN,$FP,$FN,$Sn,$Sp) if(wantarray);
  return {
    TP=>$TP,
    TN=>$TN,
    FP=>$FP,
    FN=>$FN,
    Sn=>$Sn,
    SP=>$Sp,
  };
}


sub lyveSetSnps{
  my($proj,$settings)=@_;
  my %pos;

  # Get actual SNPs from filteredMatrix
  my $matrix="$proj/msa/out.filteredMatrix.tsv";
  logmsg "Reading $matrix";
  open(SNPMATRIX,"<",$matrix) or die "ERROR: could not open $matrix for reading: $!";
  while(<SNPMATRIX>){
    next if(/^#/);
    chomp;
    my($chr,$pos,$ref,@alt)=split /\t/;
    $pos{$chr}{$pos}=$ref;
  }
  close SNPMATRIX;
  
  # Get nonsnps from snpmatrix
  my $allMatrix="$proj/msa/out.snpmatrix.tsv";
  logmsg "Reading $allMatrix";
  open(SNPMATRIX,"<",$allMatrix) or die "ERROR: could not open $allMatrix for reading: $!";
  while(<SNPMATRIX>){
    next if(/^#/);
    chomp;
    my($chr,$pos,$ref,@alt)=split /\t/;
    $pos{$chr}{$pos}//="";
  }
  close SNPMATRIX;

  return \%pos;
}

sub snvphylSnps{
  my($proj,$settings)=@_;
  my %pos;

  my $matrix=(glob("$proj/*-pseudo-positions.tsv.tabular"))[0];
  return {} if(!-e $matrix);

  logmsg "reading $matrix";
  open(SNPMATRIX,"<",$matrix) or die "ERROR: could not read $matrix: $!";
  my $header=<SNPMATRIX>; # burn the header line
  while(<SNPMATRIX>){
    chomp;
    my($chr,$pos,$status,$ref,@GT)=split /\t/;
    next if($status ne "valid");
    $pos{$chr}{$pos}||=$ref;
  }
  close SNPMATRIX;
  return \%pos;
}

# TODO: how many high-quality sites does snp-pipeline report?
sub snppipelineSnps{
  my($proj,$settings)=@_;
  my %pos;

  my $list="$proj/snplist.txt";
  open(SNPMATRIX,"<",$list) or die "ERROR: could not open $list for reading: $!";
  while(<SNPMATRIX>){
    next if(/^#/);
    chomp;
    my($chr,$pos,$count,@genomes)=split /\t/;
    $pos{$chr}{$pos}=1;
  }
  close SNPMATRIX;

  return \%pos;
}

sub realphySnps{
  my($proj,$settings)=@_;
  my %pos;

  for my $detailsFile(glob("$proj/reference/PolySeqOut_NoGenes/*details.txt")){
    open(my $fh, "<", $detailsFile) or die "ERROR: could not read $detailsFile: $!";
    my $header=<$fh>;
    while(my $line=<$fh>){
      chomp($line);
      my($strain,$contig,$gene,$refgenestartpos,$refgenepos,$refgenomepos,$orig,$poly)=split(/\t/,$line);
      $pos{$contig}{$refgenomepos}=$orig;
    }
    close $fh;
  }
  # After all the duplicate site counting, how many SNPs are there?
  my $numPositives=0;
  while(my($contig,$posHash)=each(%pos)){
    while(my($pos,$ref)=each(%$posHash)){
      $numPositives++;
    }
  }

  # Figure out the length of the genome
  my $cumulativeNumPassed=1e15;
  my $totalLength=0;
  my $sizeFile="$proj/reference/PolySeqOut_NoGenes/coreGenomeSize.txt";
  open(my $fh, "<" , $sizeFile) or die "ERROR: could not read $sizeFile $!";
  while(<$fh>){
    chomp;
    my($genome,$counts)=split(/\t/,$_);
    my($numMapped,$cumulativeNumPassedSoFar);
    ($numMapped,$cumulativeNumPassedSoFar)=(0) x 3;
    if($counts=~/(\d+)\/(\d+)\|(\d+)/){
      ($numMapped,$totalLength,$cumulativeNumPassedSoFar)=($1,$2,$3);
      $cumulativeNumPassed=$cumulativeNumPassedSoFar if($cumulativeNumPassedSoFar < $cumulativeNumPassed);
    }
  }
  close $fh;

  # Enter a number of fake sites that represent the negatives.
  my $numNegatives=$cumulativeNumPassed-$numPositives;
  makeFakePositions(\%pos,$numNegatives,$settings);

  return \%pos;
}

sub ksnp3Snps{
  my($proj,$settings)=@_;
  my %pos;

  # fixKsnpVcf.pl -ref reference.fasta < ksnp.vcf > fixed.vcf

  if(!$$settings{ref}){
    die "ERROR: need --reference for ksnp3 SNPs";
  }

  my $vcf="$proj/kSNP3.vcf";
  open(my $fh, "cat $vcf | fixKsnpVcf.pl -ref $$settings{ref} | ") or die "ERROR: could not read $vcf: $!";
  while(<$fh>){
    next if(/^#/);
    chomp;
    my($chrom,$pos,$id,$ref,$alt,$qual,$filter,$info,$format,@values)=split /\t/;
    $pos{$chrom}{$pos}=$ref;
  }
  close $fh;

  return \%pos;
}

#sub ttrSites{
#  my($proj,$settings)=@_;
#  my %pos;
#
#  my $matrix="$proj/mutsites.txt";
#  open(TTRMATRIX,"<",$matrix) or die "ERROR: could not open $matrix for reading: $!";
#  while(<TTRMATRIX>){
#    chomp;
#    my($pos)=split /\s+/;
#    $pos{$pos}=1
#  }
#  close TTRMATRIX;
#  return \%pos;
#}

sub ttrSnps{
  my($proj,$settings)=@_;
  my %pos;

  my $matrix="$proj/var_site_matrix";
  open(TTRMATRIX,"<",$matrix) or die "ERROR: could not open $matrix for reading: $!";
  while(<TTRMATRIX>){
    chomp;
    my($chr,$ref,$pos)=split /\s+/;
    $pos{$pos}=$ref;
  }
  close TTRMATRIX;

  logmsg "TODO: find how large the genome size is from the TTR directory";

  return \%pos;
}

sub makeFakePositions{
  my($pos, $numPositions, $settings)=@_;
  # Fake N positions in %pos to have no SNPs. Make them on
  # a fake contig.
  my $contigNum=1;
  my $fakeContig="FAKE_${contigNum}";
  while(defined($$pos{$fakeContig})){
    $contigNum++;
    $fakeContig="FAKE_${contigNum}";
  }
  for(my $i=1; $i<=$numPositions;$i++){
    $$pos{$fakeContig}{$i}="";
  }
}

sub flatten{
  return map{@$_} @_;
}

sub usage{
  "Compares a Lyve-SET run to a simulated dataset
  Usage: $0 --lyveset projdirectory --ttr treetoreadsdirectory
  --lyveset      proj  This is a Lyve-SET project directory
  --ttr          proj  This is the output file of a 
                       Tree To Reads run
  --snppipeline  proj  Snp-Pipeline output directory
  --RealPhy      proj  RealPhy output directory. 
                       The proj/reference subdirectory
                       will be used.
  --ksnp3        proj  kSNP3 output directory
  --snvphyl      proj  SNVPhyl project directory
  --reference    ''    Reference.fasta file. Not needed in
                       most cases.
  "
}
