#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;
use Config::Simple;
use Getopt::Long;
use File::Spec::Functions qw/rel2abs/;
use File::Basename qw/basename/;
use Bio::TreeIO;

#my $baseDir="/scicomp/home/gzu2/projects/wgsStandards/accuracyVsCoverage/manyCoverages";
my $baseDir=".";

sub logmsg{print STDERR "@_\n";}
exit(main());

sub main{
  my $settings={};
  GetOptions($settings,qw(help outprefix=s min_coverage=i max_coverage=i reps=i referencedir=s illuminaProfileDir=s tree=s numsites=i)) or die $!;
  $$settings{reps} && logmsg "WARNING: --reps is deprecated";
  $$settings{numsites}||=1;
  $$settings{outprefix}||="cov";
  $$settings{min_coverage}||=1;
  $$settings{max_coverage}||=$$settings{min_coverage};
  $$settings{referencedir}||="$baseDir/reference";
  $$settings{illuminaProfileDir}||="$baseDir/IlluminaErrorProfiles";
  $$settings{tree}||="$baseDir/lyve-set.dnd";
  die usage() if($$settings{help});
  $$settings{config}=shift(@ARGV);
  die "ERROR: need config file\n".usage() if(!$$settings{config});

  # Read the tree and pick the first name as the reference genome name.
  my $tree=Bio::TreeIO->new(-file=>$$settings{tree})->next_tree;
  my $refName=($tree->get_leaf_nodes())[0]->id;
  $refName=~s/-//g; # TTR doesn't like hyphens

  # Get configuration from cfg file
  my %config;
  Config::Simple->import_from($$settings{config},\%config);
  my $fixed=fixConfigValue(\%config);
  %config=%$fixed;


  for(my $i=$$settings{min_coverage};$i<=$$settings{max_coverage};$i++){
    my $coverage=($i*1);

    # make a simulation directory
    my $simDir="$$settings{outprefix}$i";
    $simDir=rel2abs($simDir);
    system("mkdir -pv $simDir");
    
    # copy all necessary files except config
    system("sed 's/-//g' $$settings{tree} > $simDir/tree.dnd"); # TTR doesn't like hyphens
    die if $?;
    mkdir "$simDir/illuminaProfile";
    system("cp $$settings{illuminaProfileDir}/*1.txt $simDir/illuminaProfile/1.txt");
    die if $?;
    system("cp $$settings{illuminaProfileDir}/*2.txt $simDir/illuminaProfile/2.txt");
    die if $?;
    mkdir "$simDir/reference";
    system("cp -rvL $$settings{referencedir}/reference.fasta $simDir/reference/reference.fasta");
    die if $?;
    
    # generate a custom config file
    my %newConfig=%config;
    # update some paths
    $newConfig{'default.output_dir'}="$simDir/TTR";
    $newConfig{'default.coverage'}=$coverage;
    $newConfig{'default.treefile_path'}="$simDir/tree.dnd"; 
    $newConfig{'default.base_genome_path'}="$simDir/reference/reference.fasta";
    $newConfig{'default.error_model1'}="$simDir/illuminaProfile/1.txt";
    $newConfig{'default.error_model2'}="$simDir/illuminaProfile/2.txt";
    $newConfig{'default.number_of_variable_sites'}=$$settings{numsites};
    $newConfig{'default.base_genome_name'}=$refName;
    
    open(CFG,">","$simDir/TTR.cfg") or die "ERROR: could not open $simDir/TTR.cfg for writing: $!";
    while(my($key,$value)=each(%newConfig)){
      $key=~s/^default\.//;
      if(ref($value) eq 'ARRAY'){
        $value=join(",",@$value);
      }
      print CFG "$key = $value\n";
    }
    close CFG;
  }
}

sub fixConfigValue{
  my($value)=@_;
  if(ref($value) eq 'HASH'){
    while(my($k,$v)=each(%$value)){
      $$value{$k}=fixConfigValue($v);
    }
  } elsif(ref($value) eq 'ARRAY'){
    for(my $i=0;$i<@$value;$i++){
      $$value[$i]=fixConfigValue($$value[$i]);
    }
  } else {
    $value=~s/#.*$//;
    $value=~s/^\s+|\s+$//g;
  }
  return $value;
}

sub usage{
  "Create a bunch of different treetoreads projects
  Usage: $0 treetoreads.cfg
 #--reps         1      Number of repetitions (deprecated)
  --min_coverage 1
  --max_coverage 1
  --outprefix    cov    Output directories will have the
                        format of, e.g., cov50, representing
                        the genome coverage.
  --referencedir        The directory containing reference.fasta
  --illuminaProfileDir  The directory containing something1.txt
                        and something2.txt which is the 
                        ART Illumina profile
  --tree                The newick file to model from
  --numsites     1      How many SNP sites
  "
}
