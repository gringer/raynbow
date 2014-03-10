#!/usr/bin/perl

use warnings;
use strict;

use Pod::Usage; ## uses pod documentation in usage code
use Getopt::Long qw(:config auto_version auto_help);
use File::Basename; ## for parsing file names
use IPC::Open3;
use IO::Compress::Gzip qw(gzip $GzipError);

our $VERSION = "0.1";
our $DEBUG = 0; # reduces the length of some time-consuming tasks

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
  printf(STDERR "Generating index from contig file '%s'... ",
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
  print(STDERR "done [created '$indexBase']\n");
  return($indexBase);
}

=head2 getContigLengths(baseName)

Retrieve contig lengths from the bowtie index with base I<baseName>.

=cut

sub getContigLengths{
  my ($indexBase) = @_;
  printf(STDERR "Retrieving index lengths from '%s'... ",
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
  printf(STDERR "done [found %d contigs]\n",
         scalar(keys(%{$contigLengths})));
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
  my $fragLimit = ($DEBUG?1000:10000);
  printf(STDERR "Estimating fragment size using first $fragLimit ".
         "concordant reads from '%s' and '%s':\n",
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
  printf(STDERR "[".(" " x 10)."] ");
  while(@fragLengths < $fragLimit){
    $_ = <$sout>;
    if(int(scalar(@fragLengths) * 10 / $fragLimit) > $fragProportion){
      $fragProportion = int(scalar(@fragLengths) * 10 / $fragLimit);
      printf(STDERR "\r[%s%s]",
            ("x" x $fragProportion), (" " x (10-$fragProportion)));
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
  printf(STDERR "\r[".("x" x 10)."] ");
  @fragLengths = sort({$a <=> $b} @fragLengths);
  my $medSize = @fragLengths[int(scalar(@fragLengths) / 2)];
  my $MAD = 0;
  foreach (@fragLengths){
    $MAD += abs($medSize - $_);
  }
  $MAD = $MAD / scalar(@fragLengths);
  my $child_exit_status = $? >> 8;
  printf(STDERR "done [size: %d, MAD: %d]\n", $medSize, $MAD);
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
  printf(STDERR "Removing internal concordant reads from '%s' and '%s':\n",
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
  printf(STDERR "%d reads processed [%d filtered, %d included]... ",
        $readsProcessed, $readsFiltered, $readsIncluded);
  my $outBase = "$directory/preflight_filtered";
  my $outFileNameL = "${outBase}_R1.fastq.gz";
  my $outFileNameR = "${outBase}_R2.fastq.gz";
  my $outFileL = new IO::Compress::Gzip($outFileNameL);
  my $outFileR = new IO::Compress::Gzip($outFileNameR);
  while(<$sout>){
    if(/^[^@]/){
      my @F = split(/\t/, $_, 12);
      my $filter = 0; # false
      my $RC = ($F[1] & 0x10);
      my $outReadFile = ($F[1] & 0x40) ? $outFileL : $outFileR;
      # filter out internal concordant reads
      if(($F[1] & 0x02) && ($F[6] eq "=")){
          my $leftEnd = ($F[8] > 0) ? $F[3] : $F[7];
          my $rightEnd = $leftEnd + abs($F[8]);
          my $contigEnd = $contigLengths->{$F[2]};
          $filter = (($leftEnd >= $fragMax) &&
              (($contigEnd - $rightEnd) >= $fragMax));
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
        printf(STDERR "\r%d reads processed [%d filtered, %d included]... ",
               $readsProcessed, $readsFiltered, $readsIncluded);
        if($DEBUG && $readsProcessed >= 100000){
          last;
        }
      }
    }
  }
  close($wtr);
  close($sout);
  close($serr);
  waitpid($pid, 0);
  close($outFileL);
  close($outFileR);
  printf(STDERR "\r%d reads processed [%d filtered, %d included]... ",
         $readsProcessed, $readsFiltered, $readsIncluded);
  print(STDERR "done\n");
  return(($outFileNameL, $outFileNameR, $readsIncluded));
}


=head2 edgeMap(nParallel,fragMin,fragMax,index,contigLengths,readCount,left,right)

Carries out a Bowtie2 mapping of I<left> and I<right> reads against
I<index>. Reads that map to the edges of the contigs in I<index> are
kept for later assembly.

=cut

sub edgeMap{
  my ($nParallel, $fragMin, $fragMax, $indexBase,
      $contigLengths, $numReads, $leftReads, $rightReads, $directory) = @_;
  printf(STDERR "Mapping reads from '%s' and '%s':\n",
         preDotted($leftReads), preDotted($rightReads));
  my ($wtr,$sout,$serr);
  use Symbol 'gensym'; $serr = gensym;
  my $pid = open3($wtr, $sout, $serr,
                  "bowtie2",
                  "-x",$indexBase,
                  "-1",$leftReads,
                  "-2",$rightReads,
                  "-p",$nParallel,
                  "-I",$fragMin,
                  "-X",$fragMax);

  my %attribs = ();
  my %seqs = ();
  my %contigIDs = ();
  my $nextContigID = 0;
  my $nextSeqID = 0;
  my $outFileName = "$directory/Ray_input_allcontigs.fastq.gz";
  my $outFile = new IO::Compress::Gzip($outFileName);
  my $readsProcessed = 0;
  my $readsDiscarded = 0;
  my $readsUsed = 0;
  printf(STDERR "\r%d reads processed ".
         "[%d discarded, %d used]... ",
         $readsProcessed, $readsDiscarded, $readsUsed);
  if($numReads > 0){
    printf(STDERR "%d%%", $readsProcessed * 100 / $numReads);
  }
  while(<$sout>){
    #     0    1     2     3    4     5     6     7    8   9   10    11
    # Query,Flag,R1ref,R1pos,mapQ,Cigar,R2ref,R2pos,TLen,Seq,Qual,Other
    chomp;
    if(!/^@/){
      my @F = split(/\t/, $_, 12);
      my $included = 0; # false
      my $mapped = !($F[1] & 0x04); # 0x04 -- unmapped
      my $RC = ($F[1] & 0x10); # 0x10 -- reverse complement match
      # reverse complement mapping might span the start of the contig
      # forward mapping might span the end of the contig
      my $edgeMapped = $mapped &&
        ($RC ? $F[3] <= $fragMax :
         ($contigLengths->{$F[2]} - $F[3]) <= $fragMax);
      my $edgeMappedStr = $edgeMapped ? ($RC?"S":"E") : "-";
      # printf(STDERR "%s[R%d]: %s [%s,%s,%s]\n",
      #        $F[0],($attribs{$F[0]}{"seen"})?"2":"1", $edgeMappedStr,
      #       $F[1],$F[2],$F[3]);
      if(!$contigIDs{$F[2]}){
        $contigIDs{$F[2]} = $nextContigID++;
      }
      $readsProcessed++;
      if(!$attribs{$F[0]}{"seen"}){ # this is the first end seen
        $attribs{$F[0]}{"seen"} = $edgeMappedStr;
        $attribs{$F[0]}{"ref"} = $F[2];
        $attribs{$F[0]}{"seq"} = $F[9];
        $attribs{$F[0]}{"qstr"} = $F[10];
      } else { # this is the second end seen
        if($edgeMapped){ # this read is edge mapped
          # add both reads to file associated with edge for this contig
          printf($outFile "@%s [%s]\n%s\n+\n%s\n@%s [%s]\n%s\n+\n%s\n",
                 $nextSeqID++,$contigIDs{$F[2]}.$edgeMappedStr,
                 $F[9],$F[10],
                 $nextSeqID++,$contigIDs{$F[2]}.$edgeMappedStr,
                 $attribs{$F[0]}{"seq"},$attribs{$F[0]}{"qstr"});
          $included = 1; # true
        }
        if($attribs{$F[0]}{"seen"} ne "-"){ # other read is edge mapped
          my $otherRef = $attribs{$F[0]}{"ref"};
          my $otherMStr = $attribs{$F[0]}{"seen"};
          # add both reads to file associated with edge for other contig
          printf($outFile "@%s [%s]\n%s\n+\n%s\n@%s [%s]\n%s\n+\n%s\n",
                 $nextSeqID++,$contigIDs{$otherRef}.$otherMStr,
                 $F[9],$F[10],
                 $nextSeqID++,$contigIDs{$otherRef}.$otherMStr,
                 $attribs{$F[0]}{"seq"},$attribs{$F[0]}{"qstr"});
          $included = 1; # true
        }
        delete($attribs{$F[0]});
        if($included){
          $readsUsed++;
        } else {
          $readsDiscarded++;
        }
        if(($readsProcessed % 10000) == 0){
          printf(STDERR "\r%d reads processed ".
                 "[%d discarded, %d used]... ",
                 $readsProcessed, $readsDiscarded, $readsUsed);
          if($numReads > 0){
            printf(STDERR "%d%%", $readsProcessed * 100 / $numReads);
          }
          if($DEBUG && $readsProcessed >= 100000){
            last;
          }
        }
      }
    }
  }
  close($wtr);
  close($sout);
  close($serr);
  waitpid($pid, 0);
  my $child_exit_status = $? >> 8;
  close($outFile);
  printf(STDERR "\r%d reads processed ".
         "[%d discarded, %d used]... ",
         $readsProcessed, $readsDiscarded, $readsUsed);
  print(STDERR "done\n");
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
  printf(STDERR "Using pre-defined fragment range of [%d,%d]bp\n",
         $options->{"fragMin"},$options->{"fragMax"});
}

my ($newLeft, $newRight, $newReadCount) =
  removeInternalReads($options->{"numCPUs"},
                      $options->{"fragMin"},$options->{"fragMax"},
                      $indexBase, $contigLengths,
                      $options->{"leftReads"},$options->{"rightReads"},
                      $options->{"outDir"});

edgeMap($options->{"numCPUs"},
        $options->{"fragMin"},$options->{"fragMax"},
        $indexBase, $contigLengths, $newReadCount,
        $newLeft, $newRight,
        $options->{"outDir"});

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
