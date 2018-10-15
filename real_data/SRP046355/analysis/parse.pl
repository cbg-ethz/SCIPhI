#!/usr/bin/env perl -w
# SCIPhI: Single-cell mutation identification via phylogenetic inference
#
# Copyright (C) 2018 ETH Zurich, Jochen Singer
#
# This file is part of SCIPhI.
#
# SCIPhI is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# SCIPhI is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with SCIPhI. If not, see <http://www.gnu.org/licenses/>.
#
# @author: Francesco Marass
#
# Pileup parser.
#
# 4 June 2017

use strict;
use Getopt::Std;

sub usage {
  print STDERR
"Usage: perl parse.pl [options] <pileup>

Options:
   -h         Print header.
   -m         Merge forward and reverse read information.
   -q <min>   Minimum base quality to consider {0 .. 41}.
              Note that it does not apply to indels.\n\n";
  exit 1;
}

my %opts;
getopts("hmq:", \%opts);

&usage if(scalar @ARGV != 1);

# check FASTQ quality encoding
my $bq_check = 1;
# minimum base quality + offset (default offset is 33)
my $bq_min   = 33;
# whether to merge fwd and rev counts
my $merge    = 0;
# print header?
my $header   = 0;

$merge  = 1 if(exists($opts{'m'}));
$header = 1 if(exists($opts{'h'}));
if(exists($opts{'q'})){
  if($opts{'q'} =~ /^(-?\d+)$/){
    if($1 >= 0 && $1 <= 41){
      $bq_min += $1;
    }else{
      print STDERR "Error: minimum base quality out of bounds: ",
        $1, " vs {0 .. 41}\n";
      exit 1;
    }
  }else{
    &usage;
  }
}
$bq_min = chr($bq_min);

# input handling
my $fh;
my $line;
my @flds;
# fwd, rev
my @depths;
# fwd a c g t n, rev a c g t n
my @counts;
# fwd indels, rev indels
my @indels;
#
my @seq;
my @qual;

# indices
my $i;
my $j;
my $k;
my $l;

# indel processing
my %h;
my $indels0;
my $indels1;


open($fh, $ARGV[0]) or die($!);
if($header){
  if($merge){
    print "chr\tpos\tref\tdepth\tA\tC\tG\tT\tN\tindels\n";
  }else{
    print "chr\tpos\tref\tdepth fwd\tdepth rev\tA\tC\tG\tT\tN\tindels fwd\ta\tc\tg\tt\tn\tindels rev\n";
  }
}
while($line = <$fh>){
  chomp($line);

  # zero counts
  @depths = (0) x 2;
  @counts = (0) x 10;
  @indels = ({}, {});

  @flds = split(/\s+/, $line);
  if(scalar @flds != 6){
    # samtools version 1.3.1 outputs 4 fields if depth is 0
    if(scalar @flds >= 4 && $flds[3] == 0){
      if($merge){
        print join("\t", @flds[0 .. 2],
          0, 0, 0, 0, 0, 0, ""), "\n";
      }else{
        print join("\t", @flds[0 .. 2], 0, 0,
          0, 0, 0, 0, 0, "",
          0, 0, 0, 0, 0, ""), "\n";
      }
      next;
    }else{
      print STDERR "Wrong number of fields found: ", scalar @flds, " (want 6)\n";
      close($fh);
      exit 1;
    }
  }
  # check FASTQ encoding
  if($bq_check){
    if($flds[5] =~ /[K-Z\[\\\]\^_`a-h]/){
      $bq_check = 0;
      # offset is 64, not 33
      $bq_min = chr(ord($bq_min) + 31);
    }elsif($flds[5] =~ /[0-9!"#$%&'()*+,-.\/:;<=>?\@A]/){
      # offset is 33
      $bq_check = 0;
    }
  }

  $i = 0;
  $j = 0;

  @seq  = split(//, $flds[4]);
  @qual = split(//, $flds[5]);

  while($i < scalar @seq){
    # print STDERR "(", $i, ", ", $j, ") => (", $seq[$i],", ", $qual[$j],")\n";
    if($seq[$i] eq '.'){
      $depths[0]++ if($qual[$j] ge $bq_min);
    }elsif($seq[$i] eq ','){
      $depths[1]++ if($qual[$j] ge $bq_min);
    }elsif($seq[$i] eq 'A'){
      if($qual[$j] ge $bq_min){
        $depths[0]++;
        $counts[0]++;
      }
    }elsif($seq[$i] eq 'C'){
      if($qual[$j] ge $bq_min){
        $depths[0]++;
        $counts[1]++;
      }
    }elsif($seq[$i] eq 'G'){
      if($qual[$j] ge $bq_min){
        $depths[0]++;
        $counts[2]++;
      }
    }elsif($seq[$i] eq 'T'){
      if($qual[$j] ge $bq_min){
        $depths[0]++;
        $counts[3]++;
      }
    }elsif($seq[$i] eq 'N'){
      if($qual[$j] ge $bq_min){
        $depths[0]++;
        $counts[4]++;
      }
    }elsif($seq[$i] eq 'a'){
      if($qual[$j] ge $bq_min){
        $depths[1]++;
        $counts[5]++;
      }
    }elsif($seq[$i] eq 'c'){
      if($qual[$j] ge $bq_min){
        $depths[1]++;
        $counts[6]++;
      }
    }elsif($seq[$i] eq 'g'){
      if($qual[$j] ge $bq_min){
        $depths[1]++;
        $counts[7]++;
      }
    }elsif($seq[$i] eq 't'){
      if($qual[$j] ge $bq_min){
        $depths[1]++;
        $counts[8]++;
      }
    }elsif($seq[$i] eq 'n'){
      if($qual[$j] ge $bq_min){
        $depths[1]++;
        $counts[9]++;
      }
    }elsif($seq[$i] eq '^'){
      # read start
      # the next character is the mapping quality
      # the one after is the first base
      # so do nothing here, skip next element of seq, keep qual where it is
      $i++;
      $j--;
    }elsif($seq[$i] eq '$'){
      # read end
      # ignore
      $j--;
    }elsif($seq[$i] eq '+'){
      # look for the length
      $k = $i;
      $k++ while($seq[$k + 1] =~ /\d/);
      $l = join('', @seq[($i + 1) .. $k]);
      # check orientation and save
      if($seq[$k + 1] =~ /[A-Z]/){
        $indels[0]{join('', '+', @seq[($k + 1) .. ($k + $l)])}++;
      }else{
        $indels[1]{uc(join('', '+', @seq[($k + 1) .. ($k + $l)]))}++;
      }
      # increase i, but careful with j: indels have no quality scores
      $i = $k + $l;
      $j--;
    }elsif($seq[$i] eq '-'){
      $k = $i;
      $k++ while($seq[$k + 1] =~ /\d/);
      $l = join('', @seq[($i + 1) .. $k]);
      if($seq[$k + 1] =~ /[A-Z]/){
        $indels[0]{join('', '-', @seq[($k + 1) .. ($k + $l)])}++;
      }else{
        $indels[1]{uc(join('', '-', @seq[($k + 1) .. ($k + $l)]))}++;
      }
      $i = $k + $l;
      $j--;
    }elsif($seq[$i] eq '*'){
      # do nothing
    }else{
      print STDERR "Error: unexpected sequence character (", $seq[$i], ")\n";
      close($fh);
      exit 1;
    }

    $i++;
    $j++;
  }


  ## output
  # reformat indels and output
  if($merge){
    %h = ();
    for $k (keys %{$indels[0]}){
      $h{$k} = $indels[0]{$k};
    }
    for $k (keys %{$indels[1]}){
      $h{$k} += $indels[1]{$k};
    }
    $indels0 = "";
    for $k (keys %h){
      $indels0 .= $k . 'x' . $h{$k} . ' ';
    }
    chop($indels0);

    print join("\t", @flds[0 .. 2],
      $depths[0] + $depths[1],
      $counts[0] + $counts[5], $counts[1] + $counts[6],
      $counts[2] + $counts[7], $counts[3] + $counts[8],
      $counts[4] + $counts[9], $indels0), "\n";
  }else{
    $indels0 = "";
    for $k (keys %{$indels[0]}){
      $indels0 .= $k . 'x' . $indels[0]{$k} . ' ';
    }
    chop($indels0);
    $indels1 = "";
    for $k (keys %{$indels[1]}){
      $indels1 .= $k . 'x' . $indels[1]{$k} . ' ';
    }
    chop($indels1);

    print join("\t", @flds[0 .. 2], $depths[0], $depths[1],
      @counts[0 .. 4], $indels0,
      @counts[5 .. 9], $indels1), "\n";
  }

}
close($fh);

