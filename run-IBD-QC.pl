#!/usr/bin/perl

## Code looks ate individual call rates stored in i.miss and output ID's of individuals with IBD > 0.185 
use strict;

my %imiss;
my %removed;

open IMISS, '<', $ARGV[0].".imiss"
        or die "Cannot open genotypes file (".$ARGV[0].".genome): $!\n";
print "Reading PLINK .imiss file ".$ARGV[0].".imiss\n";
while(<IMISS>){
	s/^\s+//;
    my @fields = split /\s+/, $_;
    $imiss{$fields[0]}{$fields[1]} = $fields[5];
}

open GENOME, '<', $ARGV[0].".genome"
        or die "Cannot open genotypes file (".$ARGV[0].".genome): $!\n";
open OUT, '>', "/home/thais/HelenaBIPMED/Rare/FlowFilter/ibd/fail-IBD-QC.txt";
print "Reading PLINK .genome file ".$ARGV[0].".genome\n";
while(<GENOME>){
    s/^\s+//;
    my @fields = split /\s+/, $_;
 	if($fields[9] > 0.125){ # Cut-off ibd > 0.125
 		if($imiss{$fields[0]}{$fields[1]}>$imiss{$fields[2]}{$fields[3]}){
 			unless($removed{$fields[0]}{$fields[1]}){
 				print OUT "$fields[0] $fields[1]\n";
 				$removed{$fields[0]}{$fields[1]} = 1;
 			}
 		}
 		elsif($imiss{$fields[0]}{$fields[1]}<$imiss{$fields[2]}{$fields[3]}){
 			unless($removed{$fields[2]}{$fields[3]}){
 				print OUT "$fields[2] $fields[3]\n";
 				$removed{$fields[2]}{$fields[3]} = 1;
 			}
 		}
 		else{
 			unless($removed{$fields[0]}{$fields[1]}){
 				print OUT "$fields[0] $fields[1]\n";
 				$removed{$fields[0]}{$fields[1]} = 1;
 			}
 		}
 	}
}
    
	

