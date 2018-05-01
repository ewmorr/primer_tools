#!/usr/bin/perl
#aa_alignment_barcode_candidates.pl
#Eric Morrison
#042818
#Usage: perl aa_alignment_barcode_candidates.pl [alignment fasta file] [min length to preserve sequence] [max degeneracy for primer]
#This script was designed to search protein sequence alignments for suitable sequences to use as a clade specific barcode for meta-barcoding studies. The clade of bacteria in question is an abundant and geograhically widespread genus occuring in soil for which 75 strain-level de novo genome sequences are available. The genome sequences were searched for potential orthologous genes shared across all of the strains generating a set of >800 protein sequences. The >800 protein sequences were aligned using clustal omega, and this script was then used to search the protein alignments for potential barcode regions. The criteria used for the search are: low degeneracy within a 20 bp primer region, no degeneracy in the two 3' nucleotides of potential forward or reverse primers, no gaps in barcode region, highest possible number of nucleotide combinations within the barcode region for good ability to resolve strains based on the barcode. Primer degeneracy and number of possible nucleotide combinations within the barcode are calculated based on the number of codons encoding each amino acid. For calculation of degeneracy in 3' nucleotides only the first two bases of the amino acid are considered.
#The script takes an amino acid alignment file as the first argument (tested with clustalo alignment), and minimum barcode length and maximum primer degeneracy as the second and third arguments, respectively. Suggested values for minimum barcode length is >100 amino acids which is suitable for sequencing using Illumina HiSeq 2 x 250 PE sequencing chemistry. Maximum primer degeneracy is recommended to be set intially at <100, with the lowest degeneracy preffered for potential primer sets. 
#The output is tab-separated and written to STDOUT as follows: forward primer position in the sequence, reverse primer position, forward primer degeneracy, reverse primer degeneracy, length of barcode sequence, number of barcode nucleotide combinations.

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
if(defined($ARGV[2]) == 0)
	{
	die "\nYou didn't input a maximum degeneracy for primers. Exiting...\n\n"
	}

my $in = $ARGV[0];
my $minLen = $ARGV[1];
my $maxDegen = $ARGV[2];

open(IN, "$in") || die "Can't open amino acid sequence file.\n";
chomp(my @aln = <IN>);

#rm MS linefeeds if present and process fasta
if(scalar(@aln) == 1)
	{
	$aln[0] =~ s/\r|\n|\r\n/:<:<:<:<::/g;
	@aln = split(":<:<:<:<::", $aln[0]);
	}
my $aln = join(":<:<:<:<", @aln);
$aln =~ s/\n|\r|\r\n//g;
@aln = split(">", $aln);
shift @aln;

#hash for sequences indexed by fasta header
my %aln;
my $seqLen;
my $numSeqs = 0;
#split alignments to hash
for my $seq (@aln)
	{
	my @seq = split(":<:<:<:<", $seq);
	my $head = shift @seq;
	$seq = join("", @seq);
	$aln{$head} = $seq;
	$seqLen = length($seq);
	$numSeqs++;
	}

###################################
#Get primer degeneracy by position#
###################################

my @primerN;
my @firstBaseN;
my @lastBaseN;
for(my $i = 0;$i<$seqLen-6;$i++)
	{	
	foreach my $id (keys %aln)
		{
		my @seq = split("", $aln{$id});
        	#Search with 7 aa window for primers
	        my @primerChunk = @seq[$i..$i+6];
	
	        #find n
	       	my $lastBase = pop(@primerChunk);
		my $firstBase = $primerChunk[0];
	        my $aaN = find_aa_n(@primerChunk);
		my $lastBaseN = two_bp_n($lastBase);
		
		#check for gaps
		if($aaN eq "NA" || $lastBaseN eq "NA")
			{
			$primerN[$i] = "NA";
			$lastBaseN[$i] = "NA";
			$firstBaseN[$i] = "NA";
			last;
			}
	       	
		#find 5' n
		my $firstBaseN = two_bp_n($firstBase);
	        
		$aaN *= $lastBaseN;
		$primerN[$i] += $aaN;
		$lastBaseN[$i] += $lastBaseN;
		$firstBaseN[$i] += $firstBaseN;
		}
	}

########################
#Find candidate primers#
########################

#find primers that meet degeneracy requirements
#user supplied max degeneracy, and no degeneracy in last two nt of primer (i.e. == 1)
my %forward;
my %reverse;
for(my $i = 0; $i<@primerN; $i++)
	{
	if($primerN[$i] eq "NA")
		{
		next;
		}
	if($primerN[$i]/$numSeqs <= $maxDegen)
		{
		if($lastBaseN[$i]/$numSeqs == 1)
			{
			$forward{$i} = $primerN[$i]/$numSeqs;
			}
		if($firstBaseN[$i]/$numSeqs == 1)
			{
			$reverse{$i} = $primerN[$i]/$numSeqs;
			}
		}
	}

#test potential barcode sequence for min length and gaps
foreach my $f (sort{$a <=> $b} keys %forward)
	{
	REVERSE:
	foreach my $r (sort{$a <=> $b} keys %reverse)
		{
		#test length
		if($r-($f+7) < $minLen){next REVERSE;}
		#test gaps, if none calculate degeneracy
		my $barcodeN = 0;
		my $barcodeLen = 0;
		foreach my $id (keys %aln)
			{	
			my $seq = substr($aln{$id}, $f+7, $r-($f+7));
			if($seq =~ /-/)
				{
				next REVERSE;
				}else{
				my @seq = split("",$seq);
				$barcodeN += find_aa_n(@seq);
				$barcodeLen = scalar(@seq);
				}	
			}
		#print forward and reverse primer position and degeneracy, then print barcode length and degeneracy
		print "$f\t$r\t";
		printf("%.0f", $forward{$f});
		print "\t";
		printf("%.0f", $reverse{$r});
		print "\t", $barcodeLen, "\t";
		printf("%e", $barcodeN/$numSeqs); 
		print "\n";
		}
	}


######
#SUBS#
######

#find numer of nucleotide combinations for an amino acid sequence
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
        return($chunkN);
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
if(defined($nt2_tab{$aa}) == 0)
	{
	return("NA")
	}else{
	return($nt2_tab{$aa});
	}
}


sub usage(){
	print STDERR q(
Usage: perl aa_alignment_barcode_candidates.pl [alignment fasta file] [min length to preserve sequence] [max degeneracy for primer]

This script was designed to search protein sequence alignments for suitable sequences to use as a clade specific barcode for meta-barcoding studies. The clade of bacteria in question is an abundant and geograhically widespread genus occuring in soil for which 75 strain-level de novo genome sequences are available. The genome sequences were searched for potential orthologous genes shared across all of the strains generating a set of >800 protein sequences. The >800 protein sequences were aligned using clustal omega, and this script was then used to search the protein alignments for potential barcode regions. The criteria used for the search are: low degeneracy within a 20 bp primer region, no degeneracy in the two 3' nucleotides of potential forward or reverse primers, no gaps in barcode region, highest possible number of nucleotide combinations within the barcode region for good ability to resolve strains based on the barcode. Primer degeneracy and number of possible nucleotide combinations within the barcode are calculated based on the number of codons encoding each amino acid. For calculation of degeneracy in 3' nucleotides only the first two bases of the amino acid are considered.

The script takes an amino acid alignment file as the first argument (tested with clustalo alignment), and minimum barcode length and maximum primer degeneracy as the second and third arguments, respectively. Suggested values for minimum barcode length is >100 amino acids which is suitable for sequencing using Illumina HiSeq 2 x 250 PE sequencing chemistry. Maximum primer degeneracy is recommended to be set intially at <100, with the lowest degeneracy preffered for potential primer sets. 

The output is tab-separated and written to STDOUT as follows: forward primer position in the sequence, reverse primer position, forward primer degeneracy, reverse primer degeneracy, length of barcode sequence, number of barcode nucleotide combinations. The primer positions reported are the first amino acid position of each primer.

);
exit;
}


