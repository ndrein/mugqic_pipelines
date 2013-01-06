#!/usr/env/perl

=head1 NAME

I<SAMtools>

=head1 SYNOPSIS

SAMtools->mpileup()

=head1 DESCRIPTION

B<SAMtools> is a library to manipulate and compute stats off of BAMs and to call variants

Input = file_name

Output = array

=head1 AUTHOR

=head1 DEPENDENCY

B<Pod::Usage> Usage and help output.

B<Data::Dumper> Used to debbug

=cut

package SAMtools;

# Strict Pragmas
#--------------------------
use strict;
use warnings;

#--------------------------

# Dependencies
#-----------------------

# SUB
#-----------------------
sub mpileupPaired {
  my $rH_cfg      = shift;
  my $sampleName  = shift;
  my $normalBam   = shift;
  my $tumorBam    = shift;
  my $seqName     = shift;
  my $outputDir   = shift;

  return mpileupBuilder($rH_cfg, $seqName, $outputDir, $sampleName, $normalBam, $tumorBam);
}

sub mpileup {
  my $rH_cfg      = shift;
  my $sampleName  = shift;
  my $normalBam   = shift;
  my $tumorBam    = shift;
  my $seqName     = shift;
  my $outputDir   = shift;

  return mpileupBuilder($rH_cfg, $seqName, $outputDir, $sampleName, $normalBam);
}

sub mpileupBuilder {
  my $rH_cfg      = shift;
  my $seqName     = shift;
  my $outputDir   = shift;
  my $sampleName  = shift;
  my $normalBam   = shift;
  my $tumorBam    = shift;

  my $refGenome = LoadConfig::getParam($rH_cfg, 'default', 'referenceFasta');
  my $outputBCF = $outputDir.$sampleName.'.'.$seqName.'.bcf'; 

  my $command;
  $command .= 'module load mugqic/samtools/0.1.18 ;';
  $command .= ' samtools mpileup';
  $command .= ' '.LoadConfig::getParam($rH_cfg, 'mpileup', 'mpileupExtraFlags');
  $command .= ' -f '.$refGenome;
  $command .= ' -r '.$seqName;
  $command .= ' '.$normalBam;
  if(defined($tumorBam)) {
    $command .= ' '.$tumorBam;
    $command .= ' | bcftools view -T pair -bvcg - > '.$outputBCF;
  }
  else {
    $command .= ' | bcftools view -bvcg - > '.$outputBCF;
  }

  return $command;
}

sub mergeFilterBCF {
  my $rH_cfg      = shift;
  my $sampleName  = shift;
  my $bcfDir      = shift;
  my $outputDir   = shift; 
  my $rA_seqNames = shift;

  my $refGenome = LoadConfig::getParam($rH_cfg, 'default', 'referenceFasta');
  my $outputBCF = $outputDir.$sampleName.'.merged.bcf'; 
  my $outputVCF = $outputDir.$sampleName.'.merged.flt.vcf'; 

  my $bcfInputs = "";
  for my $seqName (@$rA_seqNames) {
    my $bcfFile = $bcfDir.$sampleName.'.'.$seqName.'.bcf'; 

    $bcfInputs .= $bcfFile.' ';
  }
  
  my $command;
  $command .= 'module load mugqic/samtools/0.1.18 ;';
  $command .= ' bcftools cat';
  $command .= ' '.$bcfInputs;
  $command .= ' > '.$outputBCF;
  $command .= ' ; bcftools view '.$outputBCF;
  $command .= ' | vcfutils.pl varFilter';
  $command .= ' '.LoadConfig::getParam($rH_cfg, 'mergeFilterBCF', 'varfilterExtraFlags');
  $command .= ' > '.$outputVCF;

  return $command;
}
1;
