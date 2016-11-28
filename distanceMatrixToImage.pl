#!/usr/bin/env perl

use strict;
use warnings;
use Imager;
use Imager::Fill;
use Data::Dumper;
use Getopt::Long;
use List::Util qw/sum min max/;
use File::Basename qw/basename dirname/;

use Bio::Matrix::IO;

sub logmsg{local $0=basename $0; print STDERR "$0: @_\n";}

exit(main());

sub main{
  my $settings={};
  GetOptions($settings,qw(help sort-by-dist=s sort-by-matrix=s outfile=s px-per-square=i left-margin=i max-distance=i)) or die $!;
  $$settings{'px-per-square'}||=16;
  $$settings{'left-margin'}||=200;
  $$settings{outfile}||="tmp.bmp";

  die usage() if($$settings{help});

  my $infile=$ARGV[0];

  die usage() if(!$infile);

  drawHeatmap($infile,$settings);

  return 0;
}

sub drawHeatmap{
  my($infile,$settings)=@_;

  # how tall/wide is each square?
  my $pxPerSquare=16;

  open(DIST,"<",$infile);
  my $header=<DIST>;
  chomp($header);
  my @header=split(/\t/,$header);
  my $dot=shift(@header);
  my @distances;
  my %tableIn;

  my $maxDistance=0;
  my $userSuppliedMaxDistance = !!$$settings{'max-distance'} + 0;
  if($userSuppliedMaxDistance){
    $maxDistance=$$settings{'max-distance'};
  } 

  my $minDistance=0;
  while(<DIST>){
    chomp;
    my($queryGenome,@dist)=split /\t/;
    # Fix pesky dashes
    for(@dist){
      $_=0 if($_ eq '-' || !$_);
    }
    push(@distances,@dist);

    my %F;
    @F{@header}=@dist;
    $tableIn{$queryGenome}=\%F;

    $maxDistance=max(@dist) if(! $userSuppliedMaxDistance && max(@dist) > $maxDistance);
  }
  close DIST;

  # Sorting
  if(my $centralGenome=$$settings{'sort-by-dist'}){
    my @centralGenomeMatches = grep(/\Q$centralGenome\E/, @header);
    if(! @centralGenomeMatches){
      die "ERROR: could not find $centralGenome in matrix";
    }
    if(@centralGenomeMatches > 2){
      logmsg "WARNING: $centralGenome matches multiple genomes. Just using the first out of the following: @centralGenomeMatches";
    }
    $centralGenome=$centralGenomeMatches[0];

    @header=sort{
      $tableIn{$a}{$centralGenome} <=> $tableIn{$b}{$centralGenome}
    } @header;
  }
  if($$settings{'sort-by-matrix'}){
    open(OTHERMATRIX,"<",$$settings{'sort-by-matrix'}) or die "ERROR: could not open other matrix for sorting: $!";
    my $newHeader=<OTHERMATRIX>;
    chomp $newHeader;
    my @newHeader=split(/\t/,$newHeader);
    shift(@newHeader) if($newHeader[0] eq '.');
    # Make sure the headers match up
    for my $h(@newHeader){
      if(!grep {/\Q$h\E/} @header){
        die "ERROR: I could not find header $h in $infile";
      }
    }
    @header=@newHeader;
  }

  my $maxIntensity=255;
  #my $minIntensity=32;
  #my $intensityRange=$maxIntensity-$minIntensity+1;

  # Make a new table of colors
  my %tableColors=(header=>\@header);
  for(my $i=0;$i<@header;$i++){
    next if($header[$i] eq '.');
    for(my $j=0;$j<@header;$j++){
      next if($header[$j] eq '.');
      die "Internal error when looking at the distances between $header[$i] and $header[$j]" if(!defined($tableIn{$header[$i]}{$header[$j]}));

      my $colorIntensity=$tableIn{$header[$i]}{$header[$j]} / $maxDistance * $maxIntensity;
      $colorIntensity=$maxIntensity if($colorIntensity > $maxIntensity);
      my $primaryColorIntensity=$colorIntensity;
      my $otherColorIntensity=$colorIntensity/2;
      my $color=Imager::Color->new(red=>$primaryColorIntensity, blue=>$otherColorIntensity, green=>$otherColorIntensity);
      $tableColors{$header[$i]}{$header[$j]}=$color;
    }
  }

  # Draw the heatmap. X coordinates will be i x 10;
  # Y coordinates will be j x 10.
  my $imgWidth=(scalar(@header)+1)*$$settings{'px-per-square'} + $$settings{'left-margin'};
  my $imgHeight=(scalar(@header)+1)*$$settings{'px-per-square'};
  my $img=Imager->new(xsize=>$imgWidth,ysize=>$imgHeight,channels=>4);
  $img->box(filled=>1, color=>[255,255,255]);
  for(my $i=0;$i<@header;$i++){
    for(my $j=0;$j<@header;$j++){
      my $xmin=($i+0)*$$settings{'px-per-square'} + $$settings{'left-margin'};
      my $xmax=$xmin+$$settings{'px-per-square'} + $$settings{'left-margin'};
      my $ymin=($j+0)*$$settings{'px-per-square'};
      my $ymax=$ymin+$$settings{'px-per-square'};
      #logmsg join(", ",$xmin,$ymin,$xmax,$ymax,$tableColors{$header[$i]}{$header[$j]}->rgba);
      $img->box(color=>$tableColors{$header[$i]}{$header[$j]}, xmin=>$xmin, xmax=>$xmax, ymin=>$ymin, ymax=>$ymax, filled=>1);
      if($img->errstr){
        die "ERROR: could not draw box with rgba color "
          . join(", ",$tableColors{$header[$i]}{$header[$j]}->rgba)
          . " and coordinates $xmin, $ymin, $xmax, $ymax: "
          . $img->errstr if($img->errstr);
      }
    }
  }

  # write headers along left side
  my $font=Imager::Font->new(file=>dirname($0)."/misc/DejaVuSans.ttf") or die "Could not set font";
  for(my $i=0;$i<@header;$i++){
    my $y=($i+1)*$pxPerSquare;
    my $x=0;
    $img->string(x=>$x, y=>$y, string=>$header[$i], font=>$font, size=>$pxPerSquare, aa=>1, color=>'black');
  }

  $img->write(file=>$$settings{outfile}) or die $img->errstr;
}

sub usage{
  "$0: print a heatmap of a distance matrix
  Usage: $0 -o img.bmp heatmap.tsv
  --px-per-square  16      Sets the pixels per square on the heatmap
  --left-margin    200     Sets the pixels for the left margin, to
                           show genome labels
  --outfile        tmp.bmp The extension dictates the format.
  --sort-by-dist   ''      Sort the heatmap according to the
                           distance to a genome name
  --sort-by-matrix ''      Sort the heatmap according to
                           another distance matrix's ordering
  --max-distance   0       Set the max distance to set a standard
                           color range. Useful for normalizing
                           among a set of heatmaps. If a distance
                           in a matrix is bigger than this parameter,
                           it will be increased to match the max
                           distance found.
  "
}
