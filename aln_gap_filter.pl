#!/usr/bin/perl
#aln_gap_filter.pl
#Eric Morrison
#042618
#Usage: perl aln_gap_filter.pl [alignment fasta file] [min length to preserve sequence]
#The script takes an alignment file (tested with clustalo aa alignment, though it should work fine with nucleotides and/or other fasta format alignments). The alignments are searched by positions for gaps (-). Chunks of sequence meeting the minimum length requirement are retained and output as separate fasta files (e.g. inputname.1).

use strict;
use warnings;

if(scalar(@ARGV) == 0 || $ARGV[0] eq "-h")
	{
	&usage();
	}
if(defined($ARGV[1]) == 0)
	{
	die "\nYou didn't input a minimum length. Exiting...\n\n"
	}

my $in = $ARGV[0];
my $dir = $ARGV[1];
my $minLen = $ARGV[2];
system "mkdir $dir";

open(IN, "$in") || die "Can't open amino acid sequence file.\n";
$in =~ /.*\/(.*)/;
my $fileName = $1;

chomp(my @aln = <IN>);

if(scalar(@aln) == 1)
	{
	$aln[0] =~ s/\r|\n|\r\n/:<:<:<:<::/g;
	@aln = split(":<:<:<:<::", $aln[0]);
	}

my $aln = join(":<:<:<:<", @aln);
@aln = split(">", $aln);
shift @aln;

my %aln;
my @bare_aln;

#split alignments
for my $seq (@aln)
	{
	my @seq = split(":<:<:<:<", $seq);
	my $head = shift @seq;
	$seq = join("", @seq);
	$aln{$head} = $seq;
	push(@bare_aln, $seq);
	}

#find gaps
my @gap_pos;
for(my $i = 0; $i<length $bare_aln[0]; $i++)
	{
	for(my $u = 0; $u<@bare_aln; $u++)
		{
		if(substr($bare_aln[$u], $i, 1) eq "-")
			{
			my $not_uniq = 0;
			foreach my $check (@gap_pos)
				{
				if($check == $i)
					{
					$not_uniq = 1;
					}
				}
			if($not_uniq == 0)
				{	
				push(@gap_pos, $i);
				}
			}
		}
	}
if($gap_pos[0] != 0)
	{
	unshift(@gap_pos, 0);
	}
if($gap_pos[-1] != length($bare_aln[0])-1)
	{
	push(@gap_pos, length($bare_aln[0]));
	}
	
#cut seqs
my $fileCount = 0;
print "Gap positions in $in\n";
for(my $i = 1;$i<@gap_pos;$i++)
	{
	print $gap_pos[$i], " ";
	if($gap_pos[$i]-$gap_pos[$i-1] < $minLen)
		{
		next;
		}
	open(OUT, ">$dir/$fileName".".$fileCount") || die "Can't open output\n";
	$fileCount++;
	foreach my $key (keys %aln)
		{
		print OUT ">", $key, "\n", substr($aln{$key}, $gap_pos[$i-1], $gap_pos[$i]-$gap_pos[$i-1]), "\n";
		}
	}
print "\n";







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
