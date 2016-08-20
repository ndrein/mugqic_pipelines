#!/usr/bin/env perl
###############################################################################
# Allele ratio metrics
# Jeremy Schwartzentruber 16 June 2012
# Given a vcf file with annotated read counts as input, this function outputs
# all allele ratios for variants that pass the filtering criteria specified.
# These allele ratios can be later used to make a histogram, for example, which
# is useful in identifying problems with variant calling when an unusual
# distribution of allele ratios is seen.
###############################################################################

use strict;
use warnings;
use Getopt::Long;
sub usage($);


my ($vcfInput, $variantTypes, $minReadCount, $minAltCount, $minSNVReadRatio, $minIndelReadRatio, $minQuality, $minMapQ, $filterErrors, $prevSeenThreshold, $maxMAF, $help, $verbose);
GetOptions(	'vcf=s' => \$vcfInput,
			'variantTypes=s', \$variantTypes,
			'minReadCount=i', \$minReadCount,
			'minAltCount=i', \$minAltCount,
			'minSNVReadRatio=f', \$minSNVReadRatio,
			'minIndelReadRatio=f', \$minIndelReadRatio,
			'minQ=i', \$minQuality,
			'minMapQ=i', \$minMapQ,
			'filterErrors', \$filterErrors,
			'prevSeenThreshold=i', \$prevSeenThreshold,
			'maxMAF=f', \$maxMAF,
			'help|h' => \$help,
			'v|verbose' => \$verbose) or usage("");

$help and usage("");
$vcfInput or usage("ERROR: Missing --vcf argument.\n");

open (VCF_INPUT, "<".$vcfInput) or die ("Failed to open input file: $vcfInput");

while (<VCF_INPUT>)
{
	(/^\s*#/) and next;
	chomp();

	my ($chrm, $pos, $rsID, $ref, $alt, $qual, $filter, $vcfInfo, @restOfLine) = split(/\t/);
	my ($variation, $altCount, $readCount, $mapQ, $numPrevSamples);
	my $thgMaf = 0;
	my $evsMaf = 0;
	
	($vcfInfo =~ /VT=([^;\t]+)/)		and $variation = $1;
	($vcfInfo =~ /ALTC=([^;\t]+)/)		and $altCount = $1;
	($vcfInfo =~ /RDC=([^;\t]+)/)		and $readCount = $1;
	($vcfInfo =~ /MQ=([^;\t]+)/)		and $mapQ = $1;
	($vcfInfo =~ /PSN=([^;\t]+)/)		and $numPrevSamples = $1;
	($vcfInfo =~ /THGMAF=([^;\t]+)/)	and $thgMaf = $1;
	($vcfInfo =~ /EVSMAF=([^;\t]+)/)	and $evsMaf = $1;
	(!defined $altCount or !defined $readCount) and next;
	
	($altCount > $readCount) and die "Please check your input file: the number of alternative reads is higher than the total number of reads\n"; 

	my $isIndel = (length($ref) != length($alt));
	my $altRatio = $altCount / $readCount;
	
	if ((defined $minQuality and $qual < $minQuality) or
	    (defined $minReadCount and $readCount < $minReadCount) or
	    (defined $minAltCount and $altCount < $minAltCount) or
	    (defined $minMapQ and $mapQ < $minMapQ) or
	    (defined $prevSeenThreshold and $numPrevSamples > $prevSeenThreshold) or
	    (defined $maxMAF and ($thgMaf > $maxMAF or $evsMaf > $maxMAF)) or
	    (defined $variantTypes and $variantTypes !~ /$variation/) or
	    ($filterErrors and $filter =~ /SB|MQB|BQB|EDB/) or
	    (defined $minSNVReadRatio and !$isIndel and $altRatio < $minSNVReadRatio) or
	    (defined $minIndelReadRatio and $isIndel and $altRatio < $minIndelReadRatio))
	{
		next;
	}
	
	print $altRatio."\n";
}
close(VCF_INPUT);


###############################################################################

sub usage($) 
{
	print STDERR $_[0]."\n";
    print STDERR <<EOF;

Usage:
$0 --vcf FILE > alleleRatios.txt

Simply retrieves the alt and total read counts from the annotated VCF input file and
outputs these, one per line.

--vcf               FILE   VCF input file annotated with ALTC and RDC info fields.
--variantTypes      string Skips variants unless their type matches a type in this string
--minReadCount      int    Skips variants with read count < minReadCount
--minAltCount       int    Skips variants with alt count < minAltCount
--minSNVReadRatio   float  Skips SNV variants with read ratio < minSNVReadRatio
--minIndelReadRatio float  Skips INDEL variants with read count < minIndelReadRatio
--minQ              int    Skips variants with variant quality < minQ
--minMapQ           int    Skips variants with RMS mapping quality < minMapQ
--filterSSE                Skips variants that failed one of the filters SB,MQB,BQB,EDB
--prevSeenThreshold int    Skips variants seen in > prevSeenThreshold total samples
--maxMAF            float  Skips variants with MAF > maxMAF
EOF

    exit;
}
