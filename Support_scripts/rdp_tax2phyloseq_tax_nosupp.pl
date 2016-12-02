#!/usr/bin/perl
# rdp_tax2phyloseq_tax_nosupp.pl
use strict; use warnings;
use Getopt::Long;

# Use this script to take output from RDP Classifier and make it more manageable for phyloseq {R}.

my $usage = "\n\trdp_tax2phyloseq_tax.pl [-h -o <output PATH>] -i <input TXT>\n\n";

# defaults
my $help;
my $outDir = './';
my $inFile;

GetOptions (
	'help!' => \$help,
	'o=s'   => \$outDir,
	'i=s'   => \$inFile,
) or die $usage;

die $usage unless $help or $inFile;
if ($outDir !~ /\/$/) {$outDir = "$outDir\/"}

if ($help) {print $usage}
else {
	my ($fileName) = $inFile =~ /(.+)\.txt$/;
	my $outFile = "${outDir}$fileName.phyloseq.txt";

	open OUT, ">$outFile" or die "\n\tError: cannot create $outFile\n\n";
	print OUT "\tDomain\tPhylum\tClass\tSubclass\tOrder\tSuborder\tFamily\tGenus\n";
	open TEST, ">test_file.txt" or die "\n\tError: cannot create test_file.txt\n\n";

	open INF, "<$inFile" or die "\n\tError: cannot open $inFile\n\n";
	while (<INF>) {
		my $taxon;
		my ($domain,
			$phylum,
			$class,
			$subclass,
			$order,
			$suborder,
			$family,
			$genus
			) = ('NA') x 8;
		# print "$domain\t$domain_sup\t$phylum\t$phylum_sup\t$class\t$class_sup\t$subclass\t$subclass_sup\t";
		# print "$order\t$order_sup\t$suborder\t$suborder_sup\t$family\t$family_sup\t$genus\t$genus_sup\n";
		$_ =~ s/\"//g;
		print TEST "$_";

		($taxon)    = /^(OTU\d+)/;
		($domain)   = /\d\t([\w\ \/\-]+)\tdomain/   unless ($_ !~ /domain/);
		($phylum)   = /\d\t([\w\ \/\-]+)\tphylum/   unless ($_ !~ /phylum/);
		if ($phylum eq "NA") {$phylum = "${domain}_phylum"}
		($class)    = /\d\t([\w\ \/\-]+)\tclass/    unless ($_ !~ /class/);
		if ($class eq "NA") {$class = "${phylum}_class"}
		($subclass) = /\d\t([\w\ \/\-]+)\tsubclass/ unless ($_ !~ /subclass/);
		if ($subclass eq "NA") {$subclass = "${class}_subclass"}
		($order)    = /\d\t([\w\ \/\-]+)\torder/    unless ($_ !~ /order/);
		if ($order eq "NA") {$order = "${subclass}_order"}
		($suborder) = /\d\t([\w\ \/\-]+)\tsuborder/ unless ($_ !~ /suborder/);
		if ($suborder eq "NA") {$suborder = "${order}_suborder"}
		($family)   = /\d\t([\w\ \/\-]+)\tfamily/   unless ($_ !~ /family/);
		if ($family eq "NA") {$family = "${suborder}_family"}
		($genus)    = /\d\t([\w\ \/\-]+)\tgenus/    unless ($_ !~ /genus/);
		if ($genus eq "NA") {$genus = "${family}_genus"}

		print OUT "$taxon\t$domain\t$phylum\t$class\t$subclass\t$order\t$suborder\t$family\t$genus\n";
	}
	close INF; close OUT;
}