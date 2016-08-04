#!/usr/bin/perl
# usearch7_relabel_otus.pl
use strict; use warnings;
use Getopt::Long;

# Use this script to relabel the otus FASTA labels output by usearch version 7.0.1xxx

my $usage = "\n\tusearch7_relabel_otus.pl [-h -o -q -l <label STRING>] -i <FASTA>\n\n";

# defaults
my $help;
my $outDir = './';
my $quiet;
my $label = 'OTU';
my $inFile;

GetOptions (
	'help!'  => \$help,
	'quiet!' => \$quiet,
	'l=s'    => \$label,
	'o=s'    => \$outDir,
	'i=s'    => \$inFile,
	) or die $usage;

die $usage unless $help or $inFile;
if ($outDir !~ /\/$/) {$outDir = "$outDir\/"}

if ($help) {print $usage} 
else {
	my ($outName) = $inFile =~ /(\S+)\.fas*t*a*$/;
	my $outFile = "${outDir}${outName}.relabelled.fasta";
	open INF, "<$inFile" or die "\n\tError: cannot open $inFile\n\n";
	open OUT, ">$outFile" or die "\n\tError: cannot create $outFile\n\n";

	my $seqCount = 0;
	while (<INF>) {
		if ($_ =~ /^\>/) {
			unless ($quiet) {
				if ($seqCount % 10000 == 0) {
					$|++;
					print "\r\tseqs:\t$seqCount";
				}
			}
			my $seq = <INF>;
			chomp $seq;
			print OUT "\>${label}$seqCount\n$seq\n";
			$seqCount++;
		}
	}
	print "\r\tseqs:\t$seqCount\n\n" unless $quiet;
	close INF; close OUT;
}