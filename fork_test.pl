#!/usr/bin/perl

use warnings;
use strict;

use IO::Select;

sub indexMatch{
  my ($matchVal, $aRef) = @_;
  my $nextIndex = 0;
  my @matching = ();
  for(my $i = 0; $i < scalar(@{$aRef}); $i++){
    if(${$aRef}[$i] eq $matchVal){
      push(@matching, $i);
    }
  }
  return(@matching);
}

my $numParallel = 10;
my $progress = [(".") x $numParallel];

my $fhSelector = IO::Select->new();

$\ = $/;

my $thingsToDo = $numParallel * 2;

my $loop = 0;

my $totalThings = $thingsToDo;
my $thingsDone = 0;

printf(STDERR "\rwait step %d for processes... {%s} [%d of %d tasks done]",
       $loop++,join("",@{$progress}),$thingsDone, $totalThings);
while(($thingsToDo > 0) || (indexMatch(".",$progress) < $numParallel)){
  if($thingsToDo > 0){
    foreach my $pNum (indexMatch(".",$progress)){
      $thingsToDo--;
      ${$progress}[$pNum] = "-";
      my $pid = open(my $pipe, "-|");
      if($pid == 0){ # fork to a child process
        foreach my $letter ("A".."Z"){
          print("$pNum:$letter");
          select(undef,undef,undef, rand(1)); # sleep for up to 1 sec
        }
        close($pipe);
        exit;
      }
      $fhSelector->add($pipe);
    }
  }
  select(undef,undef,undef, 0.5); # sleep for 0.5 secs
  while(my @ready = $fhSelector->can_read(0)){
    foreach my $fh (@ready){
      my $data = <$fh>;
      if($data){
        chomp($data);
        my ($id,$val) = split(/:/, $data);
        ${$progress}[$id] = $val;
        if($val eq "Z"){
          $fhSelector->remove($fh);
          close($fh);
          ${$progress}[$id] = ".";
          $thingsDone++;
        }
      }
    }
  }
  printf(STDERR "\rwait step %d for processes... {%s} [%d of %d tasks done]",
         $loop++,join("",@{$progress}),$thingsDone, $totalThings);
}

printf(STDERR "\rwait step %d for processes... {%s} [%d of %d tasks done]\n",
       $loop++,join("",@{$progress}),$thingsDone, $totalThings);
