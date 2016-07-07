#!/usr/bin/perl
# simple_fq2fa.pl
use strict; use warnings;
use Getopt::Long;

# Use this script to convert a fastq file into a fasta file

my $usage = "\n\tUsage: simple_fq2fa.pl [options: -h -o <output directory>] -i <fastq file>\n\n";

#default options
my $help;
my $outdir = './';
my $infile;

# "global" variables
my $outext = "fa";

GetOptions (
	'o=s'   => \$outdir,
	'i=s'   => \$infile,
	'help!' => \$help,

) or die $usage;

die $usage unless(defined $infile or $help);
die "\n\tOutput directory must end in \"\/\".\n\n" unless $outdir =~ /\/$/;

if ($help) {
	print $usage;
	print "\t\t-h: this helpful help screen.\n";
	print "\t\t-o: output directory, default is current directory (\"./\").\n";
	print "\t\t-i: input fastq file.\n\n";
} else {

	open CHECK, "<$infile" or die "Error: cannot open $infile\n";
	my $line = <CHECK>;
	if ($line !~ /^\@/) {
		die "\n\tInput file must be in FASTQ format.\n\n";
	}
	close CHECK;

	open IN, "<$infile" or die "\n\tError: cannot open $infile (2nd time)\n\n";
	open OUT, ">${outdir}$infile.$outext" or die "\n\tError: cannot create ${outdir}$infile.$outext\n\n";

	while (<IN>) {
		my ($id) = $_ =~ /^@(.+)/;
		$_ = <IN>;
		my ($seq) = $_;
		$_ = <IN>;
		$_ = <IN>;
		chomp $seq;
		print OUT "\>$id\n$seq\n";
	}
	close IN; close OUT;


}
