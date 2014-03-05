#!/usr/bin/perl

use warnings;
use strict;

use Pod::Usage; ## uses pod documentation in usage code
use Getopt::Long qw(:config auto_version auto_help);
use File::Basename; ## for parsing file names
use IPC::Open3;

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

=item B<-n>

Number of iterations to run (default is to keep going until no change)

=item B<-fs>

Fragment size (default guesses from initial mapping)

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
  $indexBase =~ s/\..*?$/.index/;
  my ($wtr,$sout,$serr);
  use Symbol 'gensym'; $serr = gensym;
  my $pid = open3($wtr, $sout, $serr,
                  "bowtie2-build",$contigFile, $indexBase);
  waitpid($pid, 0);
  close($wtr);
  close($sout);
  close($serr);
  print("done");
  return($indexBase);
}



my $options =
  {
   "fragSize" => "guess",
   "iterationLimit" => -1,
   "outDir" => "out_raynbow",
};

GetOptions($options,
           'fragSize|fs=i',
           'iterationLimit|n=i',
           'contigFile|c=s',
           'leftReads|1=s',
           'rightReads|2=s',
           'outDir|o=s',
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

$\ = $/; # make print add line break

## Program meat begins here

mkdir($options->{"outDir"});

makeBT2Index($options->{"contigFile"}, $options->{"outDir"});

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
