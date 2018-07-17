#!/usr/bin/perl
#Eric Morrison
#051618

use strict;
use warnings;

#SUBS
sub usage{
    print STDERR q(
Usage:
    
The script translates a sequence of amino acids into full combination of nucleotide sequences, one per fasta entry.
    
    )
}

#sub codon_table{
#
#    my %aaToNt = (
#    "I" => ("ATT", "ATC", "ATA"),
#    "L" => ("CTT", "CTC", "CTA", "CTG", "TTA", "TTG"),
#    "V" => ("GTT", "GTC", "GTA", "GTG"),
#    "F" => ("TTT", "TTC"),
#	"M" => ("ATG"),
#	"C" => ("TGT", "TGC"),
#	"A" => ("GCT", "GCC", "GCA", "GCG"),
#	"G" => ("GGT", "GGC", "GGA", "GGG"),
#	"P" =>("CCT", "CCC", "CCA", "CCG"),
#	"T" =>("ACT", "ACC", "ACA", "ACG")
#	"S" =>("TCT", "TCC", "TCA", "TCG", "AGT", "AGC")
#	"Y" =>("TAT", "TAC")
#	"W" =>("TGG")
#	"Q" =>("CAA", "CAG")
#	"N" =>("AAT", "AAC")
#	"H" =>("CAT", "CAC")
#	"E" =>("GAA", "GAG")
#	"D" =>("GAT", "GAC")
#	"K" =>("AAA", "AAG")
#	"R" =>("CGT", "CGC", "CGA", "CGG", "AGA", "AGG")
#     )
#
#}

sub codon_to_nt_table{
	my $aa = shift @_;
	my %aaToNt = (
	"I" => "ATH",
	"L" => "YTN",
	"V" => "GTN",
	"F" => "TTY",
	"M" => "ATG",
	"C" => "TGY",
	"A" => "GCN",
	"G" => "GGN",
	"P" => "CCN",
	"T" => "ACN",
	"S" => "WSN",
	"Y" => "TAY",
	"W" => "TGG",
	"Q" => "CAR",
	"N" => "AAY",
	"H" => "CAY",
	"E" => "GAR",
	"D" => "GAY",
	"K" => "AAR",
	"R" => "MGN"
	);
	return($aaToNt{$aa});
}

sub reverse_complement{
	my $seq = $_[0];
	$seq =~ tr/ATCGRKBDYMVH/TAGCYMVHRKBD/;
	$seq = reverse($seq);
	return($seq);
}
sub convert_aa_to_nt{
	my $ntSeqs = shift @_;
	my %ntSeqs = %$ntSeqs;
	foreach my $head (keys %ntSeqs){
		$ntSeqs{$head} = build_seq($ntSeqs{$head});
	}
	return(\%ntSeqs);
}

sub build_seq{
	my @seq = split("", $_[0]);
	my $seq = "";
	foreach my $aa (@seq){
		$seq .= codon_to_nt_table($aa);
	}
	return($seq);
}
#MAIN#
{
    my $in = $ARGV[0];
	open(IN, "$in") || die "Can't open input faa file\n";
	chomp(my @faa = <IN>);
	my %faa = @faa;
	
	my $ntSeqs = convert_aa_to_nt(\%faa);
	my %ntSeqs = %$ntSeqs;
	
	foreach my $headers (keys %ntSeqs){
		if($headers =~ /r$/){
			my $rcSeq = reverse_complement($ntSeqs{$headers});
			$headers.="_rc";
			$ntSeqs{$headers} = $rcSeq;
		}
		#$ntSeqs{$headers} = substr($ntSeqs{$headers}, 0, 20);
		print $headers, "\n", $ntSeqs{$headers}, "\n";
	}
}
