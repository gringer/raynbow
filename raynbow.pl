#!/usr/bin/perl

use warnings;
use strict;

use Pod::Usage; ## uses pod documentation in usage code
use Getopt::Long qw(:config auto_version auto_help);
use File::Basename; ## for parsing file names
use IPC::Open3;
use IO::Compress::Gzip qw(gzip $GzipError);

our $VERSION = "0.01";

=head1 NAME

raynbow.pl - in-silico primer walking using Ray and Bowtie2

=head1 SYNOPSIS

./raynbow.pl -c <contigs.fasta> -1 <left.reads.fq> -2 <right.reads.fq> [options]

=head2 Basic Options

=over 2

=item B<-help>

Only display this help message

=item B<-c>

Contig file in fasta format

=item B<-1>

Left reads (can be comma-separated for multiple files)

=item B<-2>

Right reads (can be comma-separated for multiple files)

=item B<-o>

Output directory (default 'out_raynbow')

=item B<-p>

Number of parallel processing threads (default 2)

=item B<-n>

Number of iterations to run (default is to keep going until no change)

=item B<-I>

Minimum fragment size (default guesses from initial mapping)

=item B<-X>

Maximum fragment size (default guesses from initial mapping)

=back

=head1 DESCRIPTION

This uses the same concept as IMAGE
[http://www.sanger.ac.uk/resources/software/pagit/], but using Ray for
the assembler (instead of Velvet) and Bowtie2 for the mapper (instead
of SMALT).

=head1 METHODS

=cut

=head2 preDotted(string)

Replace the beginning of I<string> with dots if it has length greater
than 30 characters.

=cut

sub preDotted{
  my ($s) = @_;
  $s =~ s/^(.*?)(.{30})$/...$2/;
  return($s);
}

=head2 inPath(program)

Returns true if I<program> is in the path

=cut

sub inPath{
  my ($s) = @_;
  system('which', $s);
  return(!($?));
}

=head2 makeBT2Index(file, directory)

Generates a Bowtie2 index from I<file>, placing the index in I<directory>.

=cut

sub makeBT2Index{
  my ($contigFile, $outDir) = @_;
  printf("Generating index from contig file '%s'... ",
         preDotted($contigFile));
  my $indexBase = "$outDir/".basename($contigFile);
  $indexBase =~ s/\.[^\.]+?$/.index/;
  my ($wtr,$sout,$serr);
  use Symbol 'gensym'; $serr = gensym;
  my $pid = open3($wtr, $sout, $serr,
                  "bowtie2-build",$contigFile, $indexBase);
  waitpid($pid, 0);
  my $child_exit_status = $? >> 8;
  close($wtr);
  close($sout);
  close($serr);
  print("done [created '$indexBase']");
  return($indexBase);
}

=head2 getContigLengths(baseName)

Retrieve contig lengths from the bowtie index with base I<baseName>.

=cut

sub getContigLengths{
  my ($indexBase) = @_;
  printf("Retrieving index lengths from '%s'... ",
         preDotted($indexBase));
  my ($wtr,$sout,$serr);
  use Symbol 'gensym'; $serr = gensym;
  my $contigLengths = {};
  my $pid = open3($wtr, $sout, $serr,
                  "bowtie2-inspect","-s",$indexBase);
  while(<$sout>){
    if(/^Sequence-/){
      chomp;
      my @F = split(/\t/, $_);
      $contigLengths->{$F[1]} = $F[2];
    }
  }
  waitpid($pid, 0);
  my $child_exit_status = $? >> 8;
  close($wtr);
  close($sout);
  close($serr);
  printf("done [found %d contigs]\n", scalar(keys(%{$contigLengths})));
  return($contigLengths);
}

=head2 getFragmentSize(nParallel,index,left,right,directory)

Carries out a Bowtie2 mapping of I<left> and I<right> reads against
I<index>, using a small subset of the reads to determine mean fragment
size. Only the first 10,000 concordant reads are sampled. This number
was chosen based on data from a single run using trimmomatic-trimmed
reads from a bacterial genome.

=cut

sub getFragmentSize{
  my ($nParallel, $indexBase, $leftReads, $rightReads, $directory) = @_;
  my $fragLimit = 10000;
  printf("Estimating fragment size using first $fragLimit concordant ".
         "reads from '%s' and '%s':\n",
         preDotted($leftReads), preDotted($rightReads));
  my ($wtr,$sout,$serr);
  my @fragLengths;
  use Symbol 'gensym'; $serr = gensym;
  # run Bowtie2, piping STDOUT to $sout and STDERR to $serr
  # set fragment size range to 0~10,000, assuming sequencing < 10kbp fragments
  my $pid = open3($wtr, $sout, $serr,
                  "bowtie2",
                  "-x",$indexBase,
                  "-1",$leftReads,
                  "-2",$rightReads,
                  "-p",$nParallel,
                  "-I","0",
                  "-X","10000");
  my $fragTotal = 0;
  my $fragProportion = -1;
  $| = 1; # force autoflush
  printf("[".(" " x 50)."] ");
  while(@fragLengths < $fragLimit){
    $_ = <$sout>;
    if(int(scalar(@fragLengths) * 100 / $fragLimit) > $fragProportion){
      $fragProportion = int(scalar(@fragLengths) * 50 / $fragLimit);
      printf("\r[%s%s]",
            ("x" x $fragProportion), (" " x (50-$fragProportion)));
    }
    if(/^[^@]/){
      my @F = split(/\t/, $_, 12);
      # only record concordant reads to the same contig with positive
      # fragment size calculations
      if((($F[1] & 0x02) == 0x02) && ($F[6] eq "=") && ($F[8] > 0)){
        $fragTotal += $F[8];
        push(@fragLengths, $F[8]);
      }
    }
  }
  close($wtr);
  close($sout);
  close($serr);
  waitpid($pid, 0);
  $| = 0; # disable autoflush
  printf("\r[".("x" x 50)."] ");
  @fragLengths = sort({$a <=> $b} @fragLengths);
  my $medSize = @fragLengths[int(scalar(@fragLengths) / 2)];
  my $MAD = 0;
  foreach (@fragLengths){
    $MAD += abs($medSize - $_);
  }
  $MAD = $MAD / scalar(@fragLengths);
  my $child_exit_status = $? >> 8;
  printf("done [size: %d, MAD: %d]\n", $medSize, $MAD);
  return(($medSize - $MAD, $medSize + $MAD));
}

=head2 removeInternalReads(nParallel,fragMin,fragMax,
    index,contigLengths,left,right,directory)

Carries out a Bowtie2 mapping of I<left> and I<right> reads against
I<index>, removing fragments that sit entirely within a contig (left
position greater than one fragment length from the start, and right
position greater than one fragment length from the end).

=cut

sub removeInternalReads{
  my ($nParallel, $fragMin, $fragMax, $indexBase,
      $contigLengths, $leftReads, $rightReads, $directory) = @_;
  printf("Removing internal concordant reads from '%s' and '%s':\n",
         preDotted($leftReads), preDotted($rightReads));
  my ($wtr,$sout,$serr);
  my @fragLengths;
  use Symbol 'gensym'; $serr = gensym;
  # run Bowtie2, piping STDOUT to $sout and STDERR to $serr
  my $pid = open3($wtr, $sout, $serr,
                  "bowtie2",
                  "-x",$indexBase,
                  "-1",$leftReads,
                  "-2",$rightReads,
                  "-p",$nParallel,
                  "-I",$fragMin,
                  "-X",$fragMax);
  my $readsProcessed = 0;
  my $readsFiltered = 0;
  my $readsIncluded = 0;
  $| = 1; # force autoflush
  printf("%d reads processed [%d filtered, %d included]... ",
        $readsProcessed, $readsFiltered, $readsIncluded);
  my $outBase = "$directory/preflight_filtered";
  my $outFileL = new IO::Compress::Gzip("$outBase_R1.fastq.gz");
  my $outFileR = new IO::Compress::Gzip("$outBase_R2.fastq.gz");
  while(<$sout>){
    if(/^[^@]/){
      my @F = split(/\t/, $_, 12);
      my $filter = 0; # false
      my $RC = ($F[1] & 0x10);
      my $outReadFile = ($F[1] & 0x40) ? $outFileL : $outFileR;
      # filter out internal concordant reads
      if(($F[1] & 0x02) && ($F[6] eq "=")){
          my $leftEnd = $RC ? $F[3] : $F[7];
          my $rightEnd = $leftEnd + $F[8];
          my $contigEnd = $contigLengths->{$F[2]};
          $filter = ($leftEnd <= $fragMax) || 
              (($contigEnd - $rightEnd) <= $fragMax);
      }
      $readsProcessed++;
      if($filter){
          $readsFiltered++;
      } else {
          $readsIncluded++;
          printf($outReadFile "@%s\n%s\n+\n%s\n",
              $F[0],$F[9],$F[10]);
      }
      if(($readsProcessed % 10000) == 0){
        printf("\r%d reads processed [%d filtered, %d included]... ",
               $readsProcessed, $readsFiltered, $readsIncluded);
      }
    }
  }
  close($wtr);
  close($sout);
  close($serr);
  close($outFileL);
  close($outFileR);
  waitpid($pid, 0);
  $| = 0; # disable autoflush
  printf("\r%d reads processed [%d filtered, %d included]... ",
         $readsProcessed, $readsFiltered, $readsIncluded);
  print("done");
  return(("$outBase_R1.fastq.gz","$outBase_R1.fastq.gz"));
}


=head2 doMapping(index,left,right)

Carries out a Bowtie2 mapping of I<left> and I<right> reads against I<index>.

=cut

sub doMapping{
  my ($indexBase, $leftReads, $rightReads, $directory) = @_;
  printf("Mapping reads from '%s' and '%s'... ",
         preDotted($leftReads), preDotted($rightReads));
  my ($wtr,$sout,$serr);
  use Symbol 'gensym'; $serr = gensym;
  my $pid = open3($wtr, $sout, $serr,
                  "bowtie2","--version");
  while(<$sout>){
    printf("%s",$_);
  }
  waitpid($pid, 0);
  my $child_exit_status = $? >> 8;
  close($wtr);
  close($sout);
  close($serr);
  print("done");
  return(0);
}


my $options =
  {
   "fragMin" => "guess",
   "iterationLimit" => -1,
   "numCPUs" => 2,
   "outDir" => "out_raynbow",
};

GetOptions($options,
           'fragMin|I=i',
           'fragMax|X=i',
           'iterationLimit|n=i',
           'contigFile|c=s',
           'leftReads|1=s',
           'rightReads|2=s',
           'outDir|o=s',
           'numCPUs|p=i',
) or
  pod2usage(1);

if(-e $options->{"outDir"}){
  my $iter = 0;
  while(-e $options->{"outDir"}.".$iter"){
    $iter++;
  }
  $options->{"outDir"} .= ".$iter";
}

if(!($options->{"contigFile"}) ||
   !($options->{"leftReads"}) ||
   !($options->{"rightReads"})){
  pod2usage("Error: contig file (-c) and read files (-1,-2) must be specified");
}

if(!(-f $options->{"contigFile"})){
  pod2usage("Error: specified contig file [".$options->{"contigFile"}.
            "] does not exist");
}

if(split(/,/, $options->{"leftReads"}) !=
   split(/,/, $options->{"rightReads"})){
  pod2usage("Error: number of left and right read files do not match");
}

foreach my $dir ("left","right"){
  my @readFiles = split(/,/, $options->{$dir."Reads"});
  foreach my $readFile (@readFiles){
    if(!(-f $readFile)){
      pod2usage("Error: specified $dir read file [".$readFile.
                "] does not exist");
    }
  }
}

foreach my $prog ("bowtie2","Ray"){
  if(!inPath($prog)){
    pod2usage("Error: '$prog' cannot be found in the system path. ".
              "Please install $prog.");
  }
}

if(($options->{"fragMin"} ne "guess") && (!$options->{"fragMax"})){
  pod2usage("Error: both minimum (-I) and maximum (-X) fragment size ".
            "must be specified");
}


$\ = $/; # make print add line break

## Program meat begins here

mkdir($options->{"outDir"});

my $indexBase = makeBT2Index($options->{"contigFile"}, $options->{"outDir"});
my $contigLengths = getContigLengths($indexBase);

if($options->{"fragMin"} eq "guess"){
  ($options->{"fragMin"},
   $options->{"fragMax"}) = getFragmentSize($options->{"numCPUs"},
                                           $indexBase,
                                           $options->{"leftReads"},
                                           $options->{"rightReads"});
} else {
  printf("Using pre-defined fragment range of [%d,%d]bp\n",
         $options->{"fragMin"},$options->{"fragMax"});
}

doMapping($indexBase, $options->{"leftReads"}, $options->{"rightReads"});

=head1 AUTHOR

David Eccles (gringer) 2014 <bioinformatics@gringene.org>

=head1 LICENSE

Permission to use, copy, modify, and/or distribute this software for
any purpose with or without fee is hereby granted, provided that the
above copyright notice and this permission notice appear in all
copies.

The software is provided "as is" and the author disclaims all
warranties with regard to this software including all implied
warranties of merchantability and fitness. In no event shall the
author be liable for any special, direct, indirect, or consequential
damages or any damages whatsoever resulting from loss of use, data or
profits, whether in an action of contract, negligence or other
tortious action, arising out of or in connection with the use or
performance of this software.

=head1 AVAILABILITY

The most recent version of this code can be found on github:

https://github.com/gringer/raynbow
