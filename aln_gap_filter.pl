#!/usr/bin/perl
#aln_gap_filter.pl
#Eric Morrison
#042618
#Usage: perl aln_gap_filter.pl [alignment fasta file] [min length to preserve sequence]
#The script takes an alignment file (tested with clustalo aa alignment, though it should work fine with nucleotides and/or other fasta format alignments). The alignments are searched by positions for gaps (-). Chunks of sequence meeting the minimum length requirement are retained and output as separate fasta files (e.g. inputname.1).

use strict;
use warnings;

if(scalar(@ARGV) == 0 || $ARGV[0] eq "-h"){
	&usage();
	}
if(defined($ARGV[1]) == 0){
	die "\nYou didn't input a minimum length. Exiting...\n\n"
	}

my $in = $ARGV[0];
my $dir = $ARGV[1];
my $minLen = $ARGV[2];
#system "mkdir $dir";

open(IN, "$in") || die "Can't open amino acid sequence file.\n";
$in =~ /.*\/(.*)/;
my $fileName = $1;

chomp(my @aln = <IN>);

if(scalar(@aln) == 1){
	$aln[0] =~ s/\r|\n|\r\n/:<:<:<:<::/g;
	@aln = split(":<:<:<:<::", $aln[0]);
	}

my $aln = join(":<:<:<:<", @aln);
@aln = split(">", $aln);
shift @aln;

my %aln;
my @bare_aln;

#split alignments
for my $seq (@aln){
	my @seq = split(":<:<:<:<", $seq);
	my $head = shift @seq;
	$seq = join("", @seq);
	$aln{$head} = $seq;
	push(@bare_aln, $seq);
	}

#find gaps
my @gap_pos;
for(my $i = 0; $i<length $bare_aln[0]; $i++){
	for(my $u = 0; $u<@bare_aln; $u++){
		if(substr($bare_aln[$u], $i, 1) eq "-"){
			my $not_uniq = 0;
			foreach my $check (@gap_pos){
				if($check == $i){
					$not_uniq = 1;
					}
				}
			if($not_uniq == 0){
				push(@gap_pos, $i);
				}
			}
		}
	}

#Check for no gap alns
if(scalar(@gap_pos) == 0){
	open(OUT, ">$dir/$fileName".".full") || die "Can't open output\n";
	foreach my $key (keys %aln){
		print OUT ">", $key, "_pos_0_", length($aln{$key})-1, "\n", $aln{$key}, "\n";
		}
	print "no gaps in $fileName\n";
	exit;
	}

my $fileCount = 0;

#Check for alns with only a single gap which are not properly handled by the main loop
if(scalar(@gap_pos) == 1){
	if($gap_pos[0] > $minLen){
		open(OUT, ">$dir/$fileName".".single.$fileCount") || die "Can't open output\n";
		$fileCount++;
		foreach my $key (keys %aln){
			print OUT ">", $key, "_pos_0_", $gap_pos[0]-1, "\n", substr($aln{$key}, 0, $gap_pos[0]), "\n";
			}
		}
	if(length($bare_aln[0])-$gap_pos[0] > $minLen){
		open(OUT, ">$dir/$fileName".".single.$fileCount") || die "Can't open output\n";
		$fileCount++;
		foreach my $key (keys %aln){
			print OUT ">", $key, "pos_", $gap_pos[-1]+1, "_", length($bare_aln[0])-1, "\n", substr($aln{$key}, $gap_pos[-1]+1, length($bare_aln[0])-$gap_pos[-1]-1), "\n";
			}
		}
	print "One gap at $gap_pos[0] in $fileName\n";
	exit;
	}

#cut seqs
print "Gap positions in $fileName\n";

#need to start with zero offset if it is not a gap. If there is a single gap at the first position this should find that
if($gap_pos[0] != 0){
	if($gap_pos[0] > $minLen){
		open(OUT, ">$dir/$fileName".".$fileCount") || die "Can't open output\n";
		$fileCount++;
		foreach my $key (keys %aln){
			print OUT ">", $key, "_pos_0_", $gap_pos[0]-1, "\n", substr($aln{$key}, 0, $gap_pos[0]), "\n";
			}
		}
	}

for(my $i = 1;$i<@gap_pos;$i++){
	print $gap_pos[$i], " ";
	if($gap_pos[$i]-$gap_pos[$i-1] < $minLen){
		next;
		}
	open(OUT, ">$dir/$fileName".".$fileCount") || die "Can't open output\n";
	$fileCount++;
	foreach my $key (keys %aln){
		print OUT ">", $key, "_pos_", $gap_pos[$i-1]+1, "_", $gap_pos[$i]-1, "\n", substr($aln{$key}, $gap_pos[$i-1]+1, $gap_pos[$i]-$gap_pos[$i-1]-1), "\n";
		}
	}
print "\n";

#need to end with final offset if this is not a gap. If there is a single gap at the last position this should find that
if($gap_pos[-1] != length($bare_aln[0])-1){
	if(length($bare_aln[0])-$gap_pos[-1] > $minLen){
		open(OUT, ">$dir/$fileName".".$fileCount") || die "Can't open output\n";
		$fileCount++;
		foreach my $key (keys %aln){
			print OUT ">", $key, "pos_", $gap_pos[-1]+1, "_", length($bare_aln[0])-1, "\n", substr($aln{$key}, $gap_pos[-1]+1, length($bare_aln[0])-$gap_pos[-1]-1), "\n";
			}
		}
	}







######
#SUBS#
######

sub usage(){
	print STDERR q(
Usage: perl aln_gap_filter.pl [alignment fasta file] [min length to preserve sequence]
The script takes an alignment file (tested with clustalo aa alignment, though it should work fine with nucleotides and/or other fasta format alignments). The alignments are searched by positions for gaps (-). Chunks of sequence meeting the minimum length requirement are retained and are output to separate fasta files (e.g. inputname.1).

);
exit;
}
