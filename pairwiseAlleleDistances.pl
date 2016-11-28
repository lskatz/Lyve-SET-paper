#!/usr/bin/env perl

# open an allele file from BioNumerics to uncover pairwise distances
# When comparing the sum of percentages for 1308MDGX6-1 with the similarity
# matrix from BioNumerics, this script gave a sum of 106168.11 vs their
# sum of 106071.30.  The averages were 98.21287 vs 98.21416.
# I couldn't figure out where a rounding error made the values identical
# and so I'm saying Close Enough.
# These values were derived from a lower-left hand matrix with identity
# values.
# It was difficult to directly compare since BN encodes each genome name
# in its own identifiers.

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use File::Basename;
use Math::Round qw/nearest/;
use POSIX qw/ceil floor/;

$0=fileparse $0;

sub logmsg{print STDERR "$0: @_\n";}

exit main();

sub main{
  my $settings={};
  GetOptions($settings,qw(help)) or die $!;

  my($infile)=@ARGV;
  die usage() if(!$infile || $$settings{help});
  
  my($pairwise,$numLoci)=pairwiseDistances($infile,$settings);
  my @name=keys(%$pairwise);

  # print the matrix
  print join("\t",".",@name)."\n";
  for(my $i=0;$i<@name;$i++){
    print $name[$i];
    for(my $j=0;$j<@name;$j++){
      #print "\t".join("..",$$pairwise{$name[$i]}{$name[$j]},$name[$i],$name[$j]);
      print "\t".$$pairwise{$name[$i]}{$name[$j]};
    }
    print "\n";
  }
  
  return 0;
}

sub pairwiseDistances{
  my($infile,$settings)=@_;

  my %allele;
  my @name;

  # Read the file
  open(my $fh,$infile) or die "ERROR: could not read $infile: $!";
  my $header=<$fh>;
  while(<$fh>){
    chomp;
    my($name,@allele)=split /,/;
    $allele{$name}=\@allele;
    push(@name,$name);
  }
  close $fh;

  # Compare distances
  my %distance;
  my $numLoci=scalar(@{ $allele{$name[0]} });
  # Loop through each genome against each genome
  for(my $i=0;$i<@name;$i++){
    for(my $j=0;$j<@name;$j++){
      # Initialize the counts
      my $numPairwiseDifferences=0;
      my $numPairwiseLoci=0;
      # Loop through each locus
      for(my $k=0;$k<$numLoci;$k++){
        # Do not count masked bases
        if($allele{$name[$i]}[$k] eq '?' || $allele{$name[$j]}[$k] eq '?'){
          next;
        }
        $numPairwiseLoci++;

        # It's a difference if the alleles are different.
        if($allele{$name[$i]}[$k] ne $allele{$name[$j]}[$k]){
          $numPairwiseDifferences++;
        }
      }

      # Transform from abs count to percentage
      #$distance{$name[$i]}{$name[$j]}=(1-$numPairwiseDifferences/$numPairwiseLoci)*100;
      #$distance{$name[$i]}{$name[$j]}=nearest(0.01,$distance{$name[$i]}{$name[$j]});
      $distance{$name[$i]}{$name[$j]}=$numPairwiseDifferences;
    }
  }
  logmsg "$numLoci loci";

  return (\%distance,$numLoci) if wantarray;
  return \%distance;
}

sub usage{
  "Usage: $0 alleleMatrix.csv
  "
}
