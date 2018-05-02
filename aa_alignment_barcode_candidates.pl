#!/usr/bin/perl
#aa_alignment_barcode_candidates.pl
#Eric Morrison
#042818
#Usage: perl aa_alignment_barcode_candidates.pl [alignment fasta file] [min length to preserve sequence] [max degeneracy for primer]

use strict;
use warnings;

######
#SUBS#
######
sub usage(){
	print STDERR q(
Usage: perl aa_alignment_barcode_candidates.pl [alignment fasta file] [min length to preserve sequence] [max degeneracy for primer]

This script was designed to search protein sequence alignments for suitable sequences to use as a clade specific barcode for meta-barcoding studies. The clade of bacteria in question is an abundant and geograhically widespread genus occuring in soil for which 75 strain-level de novo genome sequences are available. The genome sequences were searched for potential orthologous genes shared across all of the strains generating a set of >800 protein sequences. The >800 protein sequences were aligned using clustal omega, and this script was then used to search the protein alignments for potential barcode regions. The criteria used for the search are: low degeneracy within a 20 bp primer region, no degeneracy in the two 3' nucleotides of potential forward or reverse primers, no gaps in barcode region, highest possible number of nucleotide combinations within the barcode region for good ability to resolve strains based on the barcode. Primer degeneracy and number of possible nucleotide combinations within the barcode are calculated based on the number of codons encoding each amino acid. For calculation of degeneracy in 3' nucleotides only the first two bases of the amino acid are considered.

The script takes an amino acid alignment file as the first argument (tested with clustalo alignment), and minimum barcode length and maximum primer degeneracy as the second and third arguments, respectively. Suggested values for minimum barcode length is >100 amino acids which is suitable for sequencing using Illumina HiSeq 2 x 250 PE sequencing chemistry. Maximum primer degeneracy is recommended to be set intially at <100, with the lowest degeneracy preffered for potential primer sets. 

The output is tab-separated and written to STDOUT as follows: filename, forward primer position in the sequence, reverse primer position, forward primer degeneracy, reverse primer degeneracy, length of barcode sequence, number of barcode nucleotide combinations. The primer positions reported are the first amino acid position of each primer.

);
exit;
}
#find numer of nucleotide combinations for an amino acid sequence
sub find_aa_n{
my $aa = $_[0];
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

if(defined($aa_tab{$aa}) == 0){
	return("NA");
	}else{
	return($aa_tab{$aa});
	}
}

#Find n of the two bp leading a codon
sub two_bp_n{
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
if(defined($nt2_tab{$aa}) == 0){
	return("NA");
	}else{
	return($nt2_tab{$aa});
	}
}

sub div_by{
my $divBy = pop(@_);
my @div = @_;
for(my $i = 0; $i < @div; $i++){
	if($div[$i] eq "NA"){
		next;
		}else{
		$div[$i] /= $divBy;
		}
	}
return(@div);	
}

sub mult_arr{
my @mult = @_;
my $mult = 1;
foreach my $i (@mult){
	$mult *= $i;
	}
return($mult);
}

sub process_fasta{
my $alnRef = $_[0];
my @aln = @$alnRef;
#rm MS linefeeds
if(scalar(@aln) == 1){
	$aln[0] =~ s/\r|\n|\r\n/:<:<:<:<::/g;
	@aln = split(":<:<:<:<::", $aln[0]);
	}
my $aln = join(":<:<:<:<", @aln);
$aln =~ s/\n|\r|\r\n//g;
@aln = split(">", $aln);
shift @aln;
#split alignments to hash
my %aln;
my $seqLen;
my $numSeqs = 0;
for my $seq (@aln){
	my @seq = split(":<:<:<:<", $seq);
	my $head = shift @seq;
	$seq = join("", @seq);
	$aln{$head} = $seq;
	$seqLen = length($seq);
	$numSeqs++;
	}

my @alnVals = (\%aln, $numSeqs, $seqLen);
return(@alnVals);
}

######
#MAIN#
######
{
if(scalar(@ARGV) == 0 || $ARGV[0] eq "-h"){
	&usage();
	}
if(defined($ARGV[1]) == 0){
	die "\nYou didn't input a minimum length. Exiting...\n\n";
	}
if(defined($ARGV[2]) == 0){
	die "\nYou didn't input a maximum degeneracy for primers. Exiting...\n\n";
	}

my $in = $ARGV[0];
my $minLen = $ARGV[1];
my $maxDegen = $ARGV[2];

open(IN, "$in") || die "Can't open amino acid sequence file.\n";
chomp(my @aln = <IN>);
$in =~ /.*\/(.*)/;
my $fileName = $1;

my @alnVals = process_fasta(\@aln);
my $alnRef = shift(@alnVals);
my %aln = %$alnRef;
my $numSeqs = $alnVals[0];
my $seqLen = $alnVals[1];

my @degen;
my @degen2;
for(my $i=0;$i<$seqLen;$i++){
	foreach my $id (keys %aln){
		my $aaN = find_aa_n(substr($aln{$id}, $i, 1));
		my $aa2N = two_bp_n(substr($aln{$id}, $i, 1));
		if($aaN eq "NA"){
			$degen[$i] = "NA";
			$degen2[$i] = "NA";
			last;
			}
		$degen[$i] += $aaN;
		$degen2[$i] += $aaN;
		}
	}
@degen = div_by(@degen, $numSeqs);
@degen2 = div_by(@degen2, $numSeqs);

for(my $i = 6; $i < @degen2 - 6 - $minLen; $i++){#need to check these
	if($degen2[$i] eq "NA"){next;}
	if($degen2[$i] == 1 and mult_arr($degen[$i-6..$i-1]) < $maxDegen){
		for(my $ii = $i + $minLen; $ii < @degen2 - 6; $ii++){
			if($degen2[$ii] eq "NA"){next;}
			if($degen2[$ii] == 1 and mult_arr($degen[$ii+1..$ii+6]) < $maxDegen){
				my $seq = join("", @degen[$i..$ii]);
				if($seq =~ /NA/){next;}
				my $fwdDegen = mult_arr($degen[$i-6..$i-1]);
				my $revDegen = mult_arr($degen[$ii+1..$ii+6]);
				my $barcodeDegen = mult_arr(@degen[$i..$ii]);
				print "$fileName\t", $i-6, "\t$ii\t";
				printf("%.0f", $fwdDegen);
				print "\t";
				printf("%.0f", $revDegen);
				print "\t", $ii-$i, "\t";
				printf("%e", $barcodeDegen); 
				print "\n";
				}
			}
		}
	}
}
