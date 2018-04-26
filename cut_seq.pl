#!/usr/bin/perl

use strict;
use warnings;

my $in = $ARGV[0];
my $searchID = $ARGV[1];#seq header
$searchID =~ s/^>//;
$searchID =~ s/"//g;
my $start = $ARGV[2];
my $end = $ARGV[3];

open(IN, "$in") || die "Can't open input\n";

chomp(my @in = <IN>);

my $fasta = join("::::::::::", @in);
my @fasta = split(">", $fasta);
shift@fasta;

my %fasta;
foreach my $seq(@fasta)
	{
	my @seq = split("::::::::::", $seq);
	my$id = shift@seq;
	my@id = split(" ", $id);
	my $sequ = join("", @seq);
	if($searchID eq $id[0] || $searchID eq "all")
		{
		print ">", $id[0], "\n", substr($sequ, $start-1, $end-$start+1), "\n"; #offset starts at 0
		}
	}
	