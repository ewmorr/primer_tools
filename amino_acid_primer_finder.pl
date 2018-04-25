#!/usr/bin/perl
#Eric Morrison
#041318
#amino_acid_primer_finder.pl
#This script searches an amino acid sequence using window size of 6 or 7 aa (18 or 21 nt) and reports the number of possible nt sequence combinations for each sequence (minimum is length of nt sequence and maximum is 4^n). The search is performed by comparing aa sequence to a table with the number of nt combinations per aa.

use strict;
use warnings;


if(scalar(@ARGV) == 0 || $ARGV[0] eq "-h")
	{
	&usage();
	}


my $inaa = $ARGV[0];
my $primerLen = $ARGV[1];


open(IN, "$inaa") || die "Can't open amino acid sequence file.\n";
if(defined($ARGV[1]) == 0)
	{
	die "\nYou didn't input a primer length. Exiting...\n\n"
	}

chomp(my @aa = <IN>);

if(scalar(@aa) == 1)
	{
	$aa[0] =~ s/\r|\n|\r\n/:<:<:<:<::/g;
	@aa = split(":<:<:<:<::", $aa[0]);
	}

my $aa = join(":<:<:<:<", @aa);
@aa = split(">", $aa);
shift @aa;
foreach my $seq (@aa)
	{
	my @seq = split(":<:<:<:<", $seq);
	my $id = shift@seq;
	$seq = join("", @seq);
	@seq = split("", $seq);
	
	my $pos = 1;
	while(@seq >= $primerLen)
		{
		my @primerChunk = @seq[0..$primerLen-1];
		my $aan = find_aa_n(@primerChunk);
		print $id, "\t", $pos, "\t", $aan, "\n";
		shift @seq;
		$pos++;
		} 
	}

	
#FUNCTIONS#
sub find_aa_n{
my @chunk = @_;
my $chunkN = 1;

my %aa_tab = (
	"I" => 3,
	"L" => 6,
	"V" => 4,
	"F" => 2,
	"M" => 1,
	"C" => 2,
	"A" => 4,
	"G" => 4,
	"P" => 4,
	"T" => 4, 
	"S" => 6,
	"Y" => 2,
	"W" => 1,
	"Q" => 2,
	"N" => 2, 
	"H" => 2,
	"E" => 2,
	"D" => 2,
	"K" => 2,
	"R" => 6);
	
	foreach my $aa (@chunk)
		{
		if(defined($aa_tab{$aa}) == 0)# in case of gaps or other non-aa characters
			{
			$chunkN = "NA";
			last;
			}else{
			$chunkN *= $aa_tab{$aa};
			}
		}
	
	return($chunkN)
	}
	
sub usage(){
	print STDERR q(
Usage: perl amino_acid_primer_finder.pl [seq.faa] [primer length]
This script takes an amino acid sequence file as input. You must also designate an integer value for length of primer to serach in number amino acids (common options are 6 and 7).

The output is a list with the sequence ID, sequence position starting at the first base of the potential primer sequence, and the number of potential combinations of nucletides for that sequence.

Prints to STDOUT. Redirect like this: cmd > file.txt
);
exit;
}


	