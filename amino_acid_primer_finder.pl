#!/usr/bin/perl
#Eric Morrison
#042718
#amino_acid_primer_finder.pl
#This script searches an amino acid sequence using window size of 6+1 amino acids and reports the number of possible nt sequence combinations for each sequence, considering on ly the first 2 bases of the final aa. The search is performed by comparing aa sequence to a table with the number of nt combinations per aa. In the case of an alignment files sequences with gaps are reported as NA degeneracy. GC content of the final aa is also reported as 0, 1, or 2, for 0, 1, or 2 3' GC.
#The output is a list with the sequence ID, sequence position starting at the first base of the potential primer sequence, and the number of potential combinations of nucletides for that sequence, and GC clamp.
#Prints to STDOUT. Redirect like this: cmd > file.txt

use strict;
use warnings;


if(scalar(@ARGV) == 0 || $ARGV[0] eq "-h")
	{
	&usage();
	}


my $inaa = $ARGV[0];

open(IN, "$inaa") || die "Can't open amino acid sequence file.\n";
chomp(my @aa = <IN>);

#process MS newlines if necessary
if(scalar(@aa) == 1)
	{
	$aa[0] =~ s/\r|\n|\r\n/:<:<:<:<::/g;
	@aa = split(":<:<:<:<::", $aa[0]);
	}
#process seqs
my $aa = join(":<:<:<:<", @aa);
@aa = split(">", $aa);
shift @aa;


#do stuff. finally...
foreach my $seq (@aa)
	{
	my @seq = split(":<:<:<:<", $seq);
	my $id = shift@seq;
	$seq = join("", @seq);
	@seq = split("", $seq);
	
	my $pos = 1;
	while(@seq >= 7)
		{
		my @primerChunk = @seq[0..6];
		
		#intitalize reporting vars
		my $aan = "NA";
		my $gcClamp = "NA";
		
		#check for gaps
		my $gap = 0;
		foreach my $base(@primerChunk)
			{
			if($base eq "-")
				{
				$gap = 1;
				}
			}
		if($gap == 1)
			{
			print $id, "\t", $pos, "\t", $aan, "\t", $gcClamp, "\n";
			shift @seq;
			$pos++;
			next;
			}
		
		#find n and gc
		my $lastBase = pop(@primerChunk);

		$aan = find_aa_n(@primerChunk);
		my @lastAA = two_bp_n_GC($lastBase);	
		$aan *= $lastAA[0];
		$gcClamp = $lastAA[1];

		print $id, "\t", $pos, "\t", $aan, "\t", $gcClamp, "\n";
		shift @seq;
		$pos++;
		} 
	}


######
#SUBS#
######

#Calculate number of nt compbinations for first 6 bp
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
		$chunkN *= $aa_tab{$aa};
		}
	
	return($chunkN)
	}

#Find n and gc of aa leading nt
sub two_bp_n_GC{

my $aa = $_[0];
my %nt2_tab = ( 	
	"I" => 1,
	"L" => 2,
	"V" => 1,
	"F" => 1,
	"M" => 1,
	"C" => 1,
	"A" => 1,
	"G" => 1,
	"P" => 1,
	"T" => 1, 
	"S" => 2,
	"Y" => 1,
	"W" => 1,
	"Q" => 1,
	"N" => 1, 
	"H" => 1,
	"E" => 1,
	"D" => 1,
	"K" => 1,
	"R" => 2
	);

my %gcClamp = ( 	
	"I" => 0,
	"L" => 0,
	"V" => 0,
	"F" => 0,
	"M" => 0,
	"C" => 1,
	"A" => 2,
	"G" => 2,
	"P" => 2,
	"T" => 1, 
	"S" => 1,
	"Y" => 0,
	"W" => 1,
	"Q" => 0,
	"N" => 0, 
	"H" => 0,
	"E" => 0,
	"D" => 0,
	"K" => 0,
	"R" => 2
	);
my @aa = $nt2_tab{$aa};
push(@aa,  $gcClamp{$aa});
return(@aa);
}

#USAGE	
sub usage(){
	print STDERR q(
Usage: perl amino_acid_primer_finder.pl [seq.faa]

This script searches an amino acid sequence fasta file using window size of 6+1 amino acids and reports the number of possible nt sequence combinations for each sequence, considering on ly the first 2 bases of the final aa. The search is performed by comparing aa sequence to a table with the number of nt combinations per aa. In the case of an alignment files sequences with gaps are reported as NA degeneracy. GC content of the final aa is also reported as 0, 1, or 2, for 0, 1, or 2 3' GC.

The output is a list with the sequence ID, sequence position starting at the first base of the potential primer sequence, the number of potential combinations of nucletides for that sequence, and GC clamp.

);
exit;
}


	
