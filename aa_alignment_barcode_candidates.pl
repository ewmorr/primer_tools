#!/usr/bin/perl
#aln_gap_filter.pl
#Eric Morrison
#042618
#Usage: perl aa_alignment_barcode_candidates.pl [alignment fasta file] [min length to preserve sequence]
#The script takes an amino acid alignment file as input (tested with clustalo aa alignment). The alignments are searched by positions for gaps (-). Chunks of sequence meeting the minimum length requirement are retained and output as separate fasta files (e.g. inputname.1).

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

open(IN, "$in") || die "Can't open amino acid sequence file.\n";
chomp(my @aln = <IN>);

#strip filename to write to output
$in =~ /.*\/(.*)\.*/;
my $fileName = $1;

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
#array for keys to retain original sorting
my @alnKeys;
#split alignments to hash
for my $seq (@aln)
	{
	my @seq = split(":<:<:<:<", $seq);
	my $head = shift @seq;
	$seq = join("", @seq);
	$aln{$head} = $seq;
	push(@alnKeys, $head);
	}

#find gaps across full alignment. A gap at any position within any sequence is retained in array
my @gap_pos;
#index by sequence position
for(my $i = 0; $i<length$aln{$alnKeys[0]}; $i++)
	{
	#index through sequences
	foreach my $id (@alnKeys)
		{
		if(substr($aln{$id}, $i, 1) eq "-")
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

#put chunks of alignments in hash for degeneracy and primer search
my %alnChunks;
#Check for no gap alns
if(scalar(@gap_pos) == 0)
	{
	open(OUT, ">$dir/$fileName.full.aln") || die "Can't open output\n";
	foreach my $id (keys @alnKeys)
		{
		#store seqeunce for further processing	
		$alnChunks{"$fileName.full\t$id"} = $aln{$id};
		#print to file
		print OUT ">", $id, "_pos_0_", length($aln{$id})-1, "\n", $aln{$id}, "\n";
		}
	exit;
	}

#file counter for output alignment chunks
my $fileCount = 0;

#Check for alns with only a single gap which are not properly handled by the main loop
if(scalar(@gap_pos) == 1)
	{
	if($gap_pos[0] > $minLen)
		{
		open(OUT, ">$dir/$fileName.single.$fileCount.aln") || die "Can't open output\n";
		$fileCount++;
		foreach my $id (@alnKeys)
			{
			#store chunk for further processing
			$alnChunks{"$fileName.full\t$id"} = substr($aln{$id}, 0, $gap_pos[0]);
			print OUT ">", $id, "_pos_0_", $gap_pos[0]-1, "\n", substr($aln{$id}, 0, $gap_pos[0]), "\n";
			}
		}
	if(length($aln{$alnKeys[0]})-$gap_pos[0] > $minLen)
		{
		open(OUT, ">$dir/$fileName.single.$fileCount.aln") || die "Can't open output\n";
		$fileCount++;
		foreach my $id (@alnKeys)
			{
			#store chunk for further processing
			$alnChunks{"$fileName.full\t$id"} = substr($aln{$id}, $gap_pos[-1]+1, length($aln{$alnKeys[0]})-$gap_pos[-1]-1);
			#print chunk to file
			print OUT ">", $id, "pos_", $gap_pos[-1]+1, "_", length($aln{$alnKeys[0]})-1, "\n", substr($aln{$id}, $gap_pos[-1]+1, length($aln{$alnKeys[0]})-$gap_pos[-1]-1), "\n";
			}
		}
	exit;
	}



#Process multiple gap alignments
#start with zero offset if it is not a gap. 
if($gap_pos[0] != 0)
	{
	if($gap_pos[0] > $minLen)
		{
		open(OUT, ">$dir/$fileName.$fileCount.aln") || die "Can't open output\n";
		$fileCount++;
		foreach my $id (@alnKeys)
			{
			$alnChunks{"$fileName.full\t$id"} = substr($aln{$id}, 0, $gap_pos[0]);
			print OUT ">", $id, "_pos_0_", $gap_pos[0]-1, "\n", substr($aln{$id}, 0, $gap_pos[0]), "\n";
			}
		}
	}
#process internal gaps
for(my $i = 1;$i<@gap_pos;$i++)
	{
	if($gap_pos[$i]-$gap_pos[$i-1] < $minLen)
		{
		next;
		}
	open(OUT, ">$dir/$fileName.$fileCount.aln") || die "Can't open output\n";
	$fileCount++;
	foreach my $id (@alnKeys)
		{
		$alnChunks{"$fileName.full\t$id"} = substr($aln{$id}, $gap_pos[$i-1]+1, $gap_pos[$i]-$gap_pos[$i-1]-1);	
		print OUT ">", $id, "_pos_", $gap_pos[$i-1]+1, "_", $gap_pos[$i]-1, "\n", substr($aln{$id}, $gap_pos[$i-1]+1, $gap_pos[$i]-$gap_pos[$i-1]-1), "\n";
		}
	}
#end with final offset if this is not a gap. 
if($gap_pos[-1] != length($aln{$alnKeys[0]})-1)
	{
	if(length($aln{$alnKeys[0]})-$gap_pos[-1] > $minLen)
		{
		open(OUT, ">$dir/$fileName.$fileCount.aln") || die "Can't open output\n";
		$fileCount++;
		foreach my $id (@alnKeys)
			{
			$alnChunks{"$fileName.full\t$id"} = substr($aln{$id}, $gap_pos[-1]+1, length($aln{$alnKeys[0]})-$gap_pos[-1]-1);	
			print OUT ">", $id, "pos_", $gap_pos[-1]+1, "_", length($aln{$alnKeys[0]})-1, "\n", substr($aln{$id}, $gap_pos[-1]+1, length($aln{$alnKeys[0]})-$gap_pos[-1]-1), "\n";
			}
		}
	}







######
#SUBS#
######

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


sub usage(){
	print STDERR q(
Usage: perl aln_gap_filter.pl [alignment fasta file] [min length to preserve sequence]
The script takes an alignment file (tested with clustalo aa alignment, though it should work fine with nucleotides and/or other fasta format alignments). The alignments are searched by positions for gaps (-). Chunks of sequence meeting the minimum length requirement are retained and are output to separate fasta files (e.g. inputname.1).

);
exit;
}
