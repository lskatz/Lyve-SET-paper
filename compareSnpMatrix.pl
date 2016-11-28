#!/usr/bin/env perl
# Uses a tree to reads directory and a Snp-Pipeline run to determine
# true positives, false negatives, false positives
#

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;

exit main();

sub main{
  my $settings={};
  GetOptions($settings,qw(help union=s intersection=s)) or die $!;
  die usage() if($$settings{help} || @ARGV < 2);

  print join("\t",qw(MATRIX1 MATRIX2 TP TN FP FN))."\n";
  compareMatrices($ARGV[0],$ARGV[1],$settings);
  
  return 0;
}

sub compareMatrices{
  my($m1,$m2,$settings)=@_;

  my $pos1=readMatrix($m1,$settings);
  my $pos2=readMatrix($m2,$settings);

  # compare
  my $c=compareSnps($pos1,$pos2,$settings);

  # Report
  print join("\t",$m1,$m2,$$c{TP}, $$c{TN}, $$c{FP}, $$c{FN})."\n";

  if($$settings{union}){
    open(UNION,">",$$settings{union}) or die "ERROR: could not open $$settings{union} for writing: $!";
    for my $chr(sort{$a cmp $b} keys(%{$$c{union}})){
      for my $pos(sort{$a <=> $b} keys(%{$$c{union}{$chr}})){
        print UNION join("\t",$chr,$pos,$$c{union}{$chr}{$pos})."\n";
      }
    }
    close UNION;
  }

  if($$settings{intersection}){
    open(INTERSECTION,">",$$settings{intersection}) or die "ERROR: could not open $$settings{intersection} for writing: $!";
    for my $chr(sort{$a cmp $b} keys(%{$$c{intersection}})){
      for my $pos(sort{$a <=> $b} keys(%{$$c{intersection}{$chr}})){
        print INTERSECTION join("\t",$chr,$pos,$$c{intersection}{$chr}{$pos})."\n";
      }
    }
    close INTERSECTION;
  }
}

sub readMatrix{
  my($matrix,$settings)=@_;
  my %pos;

  open(SNPMATRIX,"<",$matrix) or die "ERROR: could not open $matrix for reading: $!";
  while(<SNPMATRIX>){
    next if(/^#/);
    chomp;
    my($chr,$pos,$ref,@genomes)=split /\t/;
    $ref||='.';
    $pos{$chr}{$pos}=$ref;
  }
  close SNPMATRIX;

  return \%pos;
}

sub compareSnps{
  my ($pos1,$pos2,$settings)=@_;

  # Initialize counts to zero
  my($TP,$TN,$FP,$FN)=split(//,"0" x 4);

  my(%union,%intersection);
  
  # How many of the real SNPs were found?
  while(my($chr,$trueSite)=each(%$pos1)){
    while(my($truePos,$trueRef)=each(%$trueSite)){
      $$pos2{$chr}//={};
      $union{$chr}{$truePos}=$trueRef;
      if($$pos2{$chr}{$truePos}){
        $TP++; # The SNP was correctly found
        $intersection{$chr}{$truePos}=$trueRef;
      } else {
        $FN++; # The SNP was found but is not a true SNP
      }
    }
  }

  # How many SNPs were found that were not real?
  while(my($chr,$site)=each(%$pos2)){
    while(my($pos,$ref)=each(%$site)){
      $$pos1{$chr}//={};
      $union{$chr}{$pos}=$ref;
      if($$pos1{$chr}{$pos}){
        # This is a true positive and was already counted in the previous loop.
      } else {
        $FP++; # A SNP was found but is not in the real set of SNPs.
      }
    }
  }

  return{
    TP=>$TP, TN=>$TN,
    FP=>$FP, FN=>$FN,
    union=>\%union,
    intersection=>\%intersection,
  };
}

sub usage{
  "Compares a SNP matrix (columns: chrom pos ref ...) to a reference
  Usage: $0 reference.tsv query.tsv
  --union         union.tsv
  --intersection  intersection.tsv
  "
}
