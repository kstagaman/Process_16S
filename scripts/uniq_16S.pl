#!/usr/bin/perl
# uniq_16S.pl
use strict; use warnings;
use Getopt::Long;

# Use this script to collapse 16S sequence reads into unique (identical) clusters
# The sample names that the reads are from must by the only thing in the FASTA ID line

my $usage = "\n\tuniq_16S.pl [-h -q -s -u -o <out PATH> -k <INTEGER>] -i <input FASTA>\n\n";

# defaults
my $help;
my $quiet;
my $size; # whether to include the string ';size=' before the abundance (for usearch) or not (for uclust)
my $print_unmatched;
my $outDir = './';
my $inFile;
my $keep = 1; 

GetOptions (
	'help!'  => \$help,
	'quiet!' => \$quiet,
	'size!'  => \$size,
	'u!'     => \$print_unmatched,
	'o=s'    => \$outDir,
	'i=s'    => \$inFile,
	'k=i'    => \$keep,
	) or die $usage;

die unless $help or $inFile;
if ($outDir !~ /\/$/) {$outDir = "$outDir\/"}

if ($help) {
	print $usage;
	print "\t\t-h: this helpful help screen\n";
	print "\t\t-q: suppress progress output\n";
	print "\t\t-s: include the string ';size=' before the abundance for usearch, omit for uclust\n";
	print "\t\t-u: print sequences with size below -s to their own file\n";
	print "\t\t-o: output directory, default is current (./)\n";
	print "\t\t-k: keep only unique clusters with abundance equal to or greater than this number, default is 1\n";
	print "\t\t-i: input FASTA file\n";
	print "\t\tUse this script to collapse 16S sequence reads into unique (identical) clusters\n";
	print "\t\tThe sample names that the reads are from and/or a unique sequence ID (separated by -) must be the ONLY things in the FASTA headers\n\n";
}
else {
	# 'global' variables
	my $spacer;
	my $hdr_style;
	if ($size) {
		$spacer = ";size=";
		$hdr_style = "usearch_hdrs";
	} else {
		$spacer = " ";
		$hdr_style = "uclust_hdrs";
	}
	my ($outName) = $inFile =~ /(.+)\.fasta$/;
	my $gt = $keep - 1; 
	my $gtkFile = "${outDir}$outName.uniq.gt${gt}.$hdr_style.fasta";
	my $gtkMap  = "${outDir}$outName.uniq_gt${gt}_map.txt";
	my $ltkFile = "${outDir}$outName.uniq.lte${gt}.$hdr_style.fasta" if $print_unmatched;
	my $ltkMap  = "${outDir}$outName.uniq_lte${gt}_map.txt"          if $print_unmatched;
	my %uniqAbunds;
	my %abundsBySmpl;
	my @uniqs;
	my @smpls;
	my $uniqID = 0;

	open INF, "<$inFile" or die "\n\tError: cannot open $inFile\n\n";
	$|++; print "\treading sequence file ... " unless $quiet;
	while (<INF>) {
		if ($_ =~ /^\>/) {
			my ($smpl) = /^\>(\w+)\-*/;
			my $seq = <INF>;
			chomp $seq;

			$uniqAbunds{$seq}++;
			$abundsBySmpl{$seq}{$smpl}++;
		}
	}
	close INF;
	$|++; print "done\n\tSorting unique sequences ... " unless $quiet;

	@uniqs = sort {$uniqAbunds{$b} <=> $uniqAbunds{$a}} keys %uniqAbunds;
	$|++; print "done\n" unless $quiet;

	open GTK, ">$gtkFile" or die "\n\tError: cannot create $gtkFile\n\n";
	open GTM, ">$gtkMap"  or die "\n\tError: cannot create $gtkMap\n\n";
	open LTK, ">$ltkFile" or die "\n\tError: cannot create $ltkFile\n\n" if $print_unmatched;
	open LTM, ">$ltkMap"  or die "\n\tError: cannot create $ltkMap\n\n"  if $print_unmatched;

	my $uniqCount = 0;
	foreach my $uniq (@uniqs) {
		$uniqCount++;

		if ($uniqAbunds{$uniq} >= $keep) {
			print GTK ">UNQ${uniqID}${spacer}$uniqAbunds{$uniq}\n$uniq\n";
			my @smpls = sort keys %{$abundsBySmpl{$uniq}};

			foreach my $smpl (@smpls) {
				print GTM "UNQ$uniqID\t$smpl\t$abundsBySmpl{$uniq}{$smpl}\n";
			}
		}
		elsif ($print_unmatched) {
			print LTK ">UNQ${uniqID}${spacer}$uniqAbunds{$uniq}\n$uniq\n";
			my @smpls = sort keys %{$abundsBySmpl{$uniq}};

			foreach my $smpl (@smpls) {
				print LTM "UNQ$uniqID\t$smpl\t$abundsBySmpl{$uniq}{$smpl}\n";
			}
		}

		unless ($quiet) {
			if ($uniqCount % 1000 == 0) {
				$|++;
				print "\r\tUnique seqs mapped: $uniqCount";
			}
		}
		$uniqID++;
	}
	$|++; print "\r\tUnique seqs mapped: $uniqCount\n\n" unless $quiet;
	close GTK; close GTM; close LTK; close LTM;
}