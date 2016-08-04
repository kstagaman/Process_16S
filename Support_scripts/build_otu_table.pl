#!/usr/bin/perl
# build_otu_table.pl
use strict; use warnings;
use Getopt::Long;

# use this script to generate an OTU table from a UCLUST .uc output file

my $usage = "\n\tusage: build_otu_table.pl [ -h -o <output PATH> -f <uclust/uparse>] -i <UC file> -m <map TXT>\n\n";

# defaults
my $help;
my $outDir = './';
my $format = 'uclust';
my $inFile;
my $mapFile;

GetOptions (
	'help!' => \$help,
	'o=s'   => \$outDir,
	'f=s'   => \$format,
	'i=s'   => \$inFile,
	'm=s'   => \$mapFile,
	) or die $usage;

if ($outDir !~ /\/$/) {$outDir = "$outDir\/"}
die unless $help or ($inFile and $mapFile);
die "\n\t-t option must be either 'uclust' or 'uparse'\n\n" unless ($format eq 'uclust' or $format eq 'uparse');

if ($help) {
	help_txt();
}
else {
	# global variables
	my ($outName) = $inFile =~ /\/*([\w\.]+)\.u[cp]$/;
	# my %seeds_by_uniqIDs;
	# my @seeds;
	# my %seqIDs_by_smpl;
	# my %abunds_by_seqID;
	# my %clusterIDs_by_seed;
	my %smplAbunds_by_uniqID;
	my %uniqIDs_by_clusterID;
	my %smplAbunds_by_clusterID;
	my @clusterIDs;
	my %uniqSmpls;
	my @smpls;
	my %uniqUniqIDs;
	my @uniqIDs;

	open MAP, "<$mapFile" or die "\n\tError: cannot open $mapFile\n\n";
	while (<MAP>) {
		my ($uniqID)    = /^(UNQ\d+)\t/;
		# print "MAP: \$uniqID = $uniqID\n";
		my ($smpl)      = /^$uniqID\t(\S+)\t/;
		# print "MAP: \$smpl = $smpl\n";
		my ($smplAbund) = /\t(\d+)$/;
		# print "MAP: \$smplAbund = $smplAbund\n";

		$smplAbunds_by_uniqID{$uniqID}{$smpl} = $smplAbund;
		# print "MAP: \$smplAbunds_by_uniqID{$uniqID}{$smpl} = $smplAbunds_by_uniqID{$uniqID}{$smpl}\n";
		$uniqSmpls{$smpl} = 0;
		# $uniqUniqIDs{$uniqID} = 0;
	}
	@smpls   = sort keys %uniqSmpls;
	# print "after MAP: @smpls\n";
	# @uniqIDs = sort keys %uniqUniqIDs;
	# print "after MAP: @uniqIDs\n";

	open INF, "<$inFile" or die "\n\tError: cannot open $inFile\n\n";
	if ($format eq 'uclust') {
		while (<INF>) {
			if ($_ =~ /^[SH]/) {  ### May or may not want to include library (L) lines
				# my ($type)       = /^([SH])\t/;  ### commented out because it went unused
				# print "INF: \$type = $type\n";
				my ($clusterID)  = /^[SH]\t(\d+)\t/;
				# print "INF: \$clusterID = $clusterID\t";
				my ($queryLabel) = /\t(UNQ\d+)\S*\t\S+$/;
				# print "INF: \$queryLabel = $queryLabel\n";

				push @{$uniqIDs_by_clusterID{$clusterID}}, $queryLabel;
				# print "\@{\$uniqIDs_by_clusterID{$clusterID}} = @{$uniqIDs_by_clusterID{$clusterID}}\n";
			}
		}
	}
	else {
		while (<INF>) {
			if ($_ =~ /\t(match|otu)\t/) { ### May or may not want to include library (???) lines
				my ($queryLabel) = /^(UNQ\d+)\;/;
				# my ($type)       = /\t(match|otu)\t/; ### commented out because it went unused
				my ($clusterID)  = /\tOTU(\d+)$/;
				push @{$uniqIDs_by_clusterID{$clusterID}}, $queryLabel;
			}
		}
	}
	close INF;

	@clusterIDs = sort {$a <=> $b} keys %uniqIDs_by_clusterID;
	open TEST, ">test.log" or die "\n\tError: cannot create test.log\n\n";
	print TEST "smpl\totu\tabund\n";
	# foreach my $clusterID (@clusterIDs) {
	# 	print TEST "$clusterID\n";
	# }
	# my @assignedUniqIDs;

	foreach my $clusterID (@clusterIDs) {
		my @uniqIDs = @{$uniqIDs_by_clusterID{$clusterID}};
		# print "\$uniqIDs_by_clusterID{$clusterID} = @uniqIDs\n";
		# push @assignedUniqIDs, @uniqIDs;
		foreach my $uniqID (@uniqIDs) {
			foreach my $smpl (@smpls) {
				# print "\$clusterID=$clusterID\t\$uniqID=$uniqID\t\$smpl=$smpl\n";
				$smplAbunds_by_uniqID{$uniqID}{$smpl} = 0 unless ($smplAbunds_by_uniqID{$uniqID}{$smpl});
				$smplAbunds_by_clusterID{$clusterID}{$smpl} += $smplAbunds_by_uniqID{$uniqID}{$smpl};
				# print "\$smplAbunds_by_uniqID{$uniqID}{$smpl} = $smplAbunds_by_uniqID{$uniqID}{$smpl}\n";
				# print "\$smplAbunds_by_clusterID{$clusterID}{$smpl} = $smplAbunds_by_clusterID{$clusterID}{$smpl}\n";
			}
		}
	}
	# @assignedUniqIDs = sort @assignedUniqIDs;
	# print "@assignedUniqIDs\n";

	my @taxa;
	foreach my $clusterID (@clusterIDs) {
		my $taxon = "OTU$clusterID";
		push @taxa, $taxon;
	}


	open OUT, ">${outDir}${outName}.otu_tbl.txt" or die "\n\tError: cannot create ${outDir}${outName}.otu_tbl.txt\n\n";
	my $columnHeader = join "\t", @taxa;
	print OUT "\t$columnHeader\n";

	# my %smplSums;
	foreach my $smpl (@smpls) {
		print OUT "$smpl";

		foreach my $clusterID (@clusterIDs) {
			print OUT "\t$smplAbunds_by_clusterID{$clusterID}{$smpl}";
			print TEST "$smpl\tX$clusterID\t$smplAbunds_by_clusterID{$clusterID}{$smpl}\n";
			# $smplSums{$smpl} += $smplAbunds_by_clusterID{$clusterID}{$smpl};
		}
		print OUT "\n";

	}
	close OUT;

	# foreach my $smpl (@smpls) {
	# 	print TEST "\$smplSums{$smpl} = $smplSums{$smpl}\n";
	# }
}

sub help_txt {
	print $usage;
}