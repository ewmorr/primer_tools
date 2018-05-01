#!/usr/bin/perl
#aa_alignment_barcode_candidates.pl
#Eric Morrison
#042818
#Usage: perl aa_alignment_barcode_candidates.pl [alignment fasta file] [min length to preserve sequence] [max degeneracy for primer]
#The script takes an amino acid alignment file as input (tested with clustalo aa alignment). The alignments are searched by positions for gaps (-). Chunks of sequence meeting the minimum length requirement are retained and output as separate fasta files (e.g. inputname.1).

use strict;
use warnings;

if(scalar(@ARGV) == 0 || $ARGV[0] eq "-h")
	{
	usage();
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

#rm MS linefeeds if present
if(scalar(@aln) == 1)
	{
	$aln[0] =~ s/\r|\n|\r\n/:<:<:<:<::/g;
	@aln = split(":<:<:<:<::", $aln[0]);
	}
#Process fasta
my $aln = join(":<:<:<:<", @aln);
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
                	$primerN[$i] = "NA";
			$lastBaseN[$i] = "NA";
			$firstBaseN[$i] = "NA";
                        last;
                        }

	        #find n
	       	my $lastBase = pop(@primerChunk);
		my $firstBase = @primerChunk[0];
	        my $aaN = find_aa_n(@primerChunk);
	        my $lastBaseN = two_bp_n($lastBase);
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
for(my $i = 0;$i<@primerN;$i++)
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
	foreach my $r (sort{$a <=> $b} keys %reverse))
		{
		#test length
		if($r-($f+7) < $minLen){next REVERSE;}
		#test gaps, if none calculate degeneracy
		my $barcodeN = 0;
		foreach my $id (keys %aln)
			{	
			my $seq = substr($aln{$id}, $f+7, $r-($f+7);
			if($seq =~ /-/)
				{
				next REVERSE;
				}else{
				my @seq = split("",$seq);
				$barcodeN += find_aa_n(@seq);
				}	
			}
		#print forward and reverse primer position and degeneracy, then print barcode length and degeneracy
		print "$f\t$r\t$forward{$f}\t$reverse{$r}\t", length($seq), $barcodeN/$numSeqs, "\n";
		}
	}


######
#SUBS#
######

sub find_aa_n(){
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
#Find n and gc of the two bp leading a codon
sub two_bp_n(){

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

return($nt2_tab{$aa});
}


sub usage(){
	print STDERR q(
Usage: perl aln_gap_filter.pl [alignment fasta file] [min length to preserve sequence]
The script takes an alignment file (tested with clustalo aa alignment, though it should work fine with nucleotides and/or other fasta format alignments). The alignments are searched by positions for gaps (-). Chunks of sequence meeting the minimum length requirement are retained and are output to separate fasta files (e.g. inputname.1).

);
exit;
}


