#!/usr/bin/perl
# join_multi-lined_fasta.pl
use strict; use warnings;

# some program output fasta files where the sequences are broken up into lines instead of being on
# 1 line.  This makes it more difficult for running scripts.  This script puts all seqs on 1 line.

my $usage = "\n\tjoin_multi-lined_fasta.pl <FASTA file>\n\n";
die $usage unless @ARGV == 1;

open IN, "<$ARGV[0]" or die "\n\tError: cannot open $ARGV[0]\n\n";
open OUT, ">$ARGV[0].joined" or die "\n\tError: cannot create $ARGV[0].joined\n\n"; 

my $first_line = <IN>;
print OUT "$first_line";

while (<IN>) {

	if ($_ !~ /^\>/) {
		my ($seq_seg) = /^([acgtnyrw\-]+)$/i;
		$seq_seg = uc $seq_seg;
		print OUT "$seq_seg";
	} else {
		my ($id) = /^\>(.+)/;
		print OUT "\n\>$id\n";
	} 
}
close IN; close OUT;