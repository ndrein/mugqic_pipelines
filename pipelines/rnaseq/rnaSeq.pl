#!/usr/bin/perl

=head1 NAME

I<rnaSeq>

=head1 SYNOPSIS

rnaSeq.pl

=head1 DESCRIPTION

B<rnaSeq> Is the main RNAseq pipeline.

=head1 AUTHOR

B<Mathieu Bourgey> - I<mbourgey@genomequebec.com>

=head1 DEPENDENCY

B<Pod::Usage> Usage and help output.

B<Data::Dumper> Used to debbug

B<Config::Simple> Used to parse config file

B<File::Basename> path parsing

B<Cwd> path parsing

=cut

# Strict Pragmas
#---------------------
use strict;
use warnings;
#---------------------

BEGIN{
    #Makesure we can find the GetConfig::LoadModules module relative to this script install
    use File::Basename;
    use Cwd 'abs_path';
    my ( undef, $mod_path, undef ) = fileparse( abs_path(__FILE__) );
    unshift @INC, $mod_path."lib";
}


# Dependencies
#--------------------
use Getopt::Std;
use Cwd qw/ abs_path /;

use LoadConfig;
use Picard;
use SampleSheet;
use SAMtools;
use SequenceDictionaryParser;
use SubmitToCluster;
use TophatBowtie;
use Trimmomatic;
use Metrics;
use Cufflinks;
use Wiggle;
use HtseqCount;
use DiffExpression;
use GqSeqUtils;

#--------------------


# SUB
#--------------------

my @steps;
push(@steps, {'name' => 'trimming' , 'stepLoop' => 'sample' , 'output' => 'reads'});
push(@steps, {'name' => 'trimMetrics' , 'stepLoop' => 'group' , 'output' => 'metrics'});
push(@steps, {'name' => 'aligning' , 'stepLoop' => 'sample' , 'output' => 'alignment'});
push(@steps, {'name' => 'merging' , 'stepLoop' => 'sample' , 'output' => 'alignment'});
push(@steps, {'name' => 'alignMetrics' , 'stepLoop' => 'group' , 'output' => 'metrics'});
#push(@steps, {'name' => 'mutation' , 'stepLoop' => 'sample' , 'output' => 'mpileup'});
push(@steps, {'name' => 'wiggle' , 'stepLoop' => 'sample' , 'output' => 'tracks'});
push(@steps, {'name' => 'rawCounts' , 'stepLoop' => 'sample' , 'output' => 'raw_counts'});
push(@steps, {'name' => 'rawCountsMetrics' , 'stepLoop' => 'group' , 'output' => 'metrics'});
push(@steps, {'name' => 'fpkm' , 'stepLoop' => 'sample' , 'output' => 'fpkm'});
push(@steps, {'name' => 'exploratory' , 'stepLoop' => 'group' , 'output' => 'exploratory'});
push(@steps, {'name' => 'cuffdiff' , 'stepLoop' => 'group' , 'output' => 'DGE'});
#push(@steps, {'name' => 'dgeMetrics' , 'stepLoop' => 'group' , 'output' => 'metrics'});
push(@steps, {'name' => 'dge' , 'stepLoop' => 'group' , 'output' => 'DGE'});
push(@steps, {'name' => 'goseq' , 'stepLoop' => 'group' , 'output' => 'DGE'});
push(@steps, {'name' => 'deliverable' , 'stepLoop' => 'group' , 'output' => 'deliverable'});



my %globalDep;
for my $stepName (@steps) { 
	$globalDep{$stepName -> {'name'} } ={};
}


# Global scope variables
my $designFilePath;
my $configFile;
my $workDir;
my $readSetSheet;


&main();

sub printUsage {
  print "\nUsage: perl ".$0." \n";
  print "\t-c  config file\n";
  print "\t-s  start step, inclusive\n";
  print "\t-e  end step, inclusive\n";
  print "\t-n  nanuq sample sheet\n";
  print "\t-d  design file\n";
  print "\t-w  work directory\n";
  print "\n";
  print "Steps:\n";
  for(my $idx = 0; $idx < @steps; $idx++) {
    print "".($idx + 1) . '- ' . $steps[$idx]->{'name'} . "\n";
  }
  print "\n";
}

sub main {
	my %opts;
	getopts('c:s:e:n:d:w:', \%opts);
	
	if (!defined($opts{'c'}) || !defined($opts{'s'}) || !defined($opts{'e'}) || !defined($opts{'n'}) || !defined($opts{'d'}) || !defined($opts{'w'} ) ) {
		printUsage();
		exit(1);
	}
	
	my %jobIdVarPrefix;
	my %cfg = LoadConfig->readConfigFile($opts{'c'});
	my $rHoAoH_sampleInfo = SampleSheet::parseSampleSheetAsHash($opts{'n'});
	my $rAoH_seqDictionary = SequenceDictionaryParser::readDictFile(\%cfg);
	$designFilePath = abs_path($opts{'d'});
	##get design groups
	my $rHoAoA_designGroup = Cufflinks::getDesign(\%cfg,$designFilePath);
	$workDir = abs_path($opts{'w'});
	$configFile =  abs_path($opts{'c'});
	$readSetSheet =  abs_path($opts{'n'}); 

	
	#generate sample jobIdprefix
	my $cpt = 1;
	
	for my $sampleName (keys %{$rHoAoH_sampleInfo}) {
		my $cpt2=1;
		$jobIdVarPrefix{$sampleName} = $cpt;
		my $rAoH_sampleLanes = $rHoAoH_sampleInfo->{$sampleName};
		for my $rH_laneInfo (@$rAoH_sampleLanes) {
			$jobIdVarPrefix{$sampleName .'.' .$rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'}} = $cpt.'_'.$cpt2;
			$cpt2++;
		}
		$cpt++;
	}
	#generate design jobIdprefix
	for my $designName (keys %{$rHoAoA_designGroup}) {
		$jobIdVarPrefix{$designName} = $cpt;
		$cpt++;
	}
	
	SubmitToCluster::initPipeline($workDir);

	my $latestBam;

	
	for(my $current = $opts{'s'}-1; $current <= ($opts{'e'}-1); $current++) {
		my $fname = $steps[$current]->{'name'};
		my $loopType = $steps[$current]->{'stepLoop'};
		my $subref = \&$fname;
		if ($loopType eq 'sample') {
			for my $sampleName (keys %{$rHoAoH_sampleInfo}) {
			my $rAoH_sampleLanes = $rHoAoH_sampleInfo->{$sampleName};
			# Tests for the first step in the list. Used for dependencies.
			my $jobIdVar = &$subref($current != ($opts{'s'}-1), \%cfg, $sampleName, $rAoH_sampleLanes, $rAoH_seqDictionary, \%jobIdVarPrefix);
			if (defined($jobIdVar)) {
				$globalDep{$fname}->{$sampleName} = $jobIdVar;
			}
			}
		}
		else {
			# Tests for the first step in the list. Used for dependencies.
			my $jobIdVar = &$subref($current != ($opts{'s'}-1), \%cfg, $rHoAoH_sampleInfo, $rHoAoA_designGroup, $rAoH_seqDictionary, \%jobIdVarPrefix);
			if (defined($jobIdVar)) {
				$globalDep{$fname}->{$fname} = $jobIdVar;
			}
		}
	}  
}

sub trimming {
	my $depends = shift;
	my $rH_cfg = shift;
	my $sampleName = shift;
	my $rAoH_sampleLanes  = shift;
	my $rAoH_seqDictionary = shift;
	my $rH_jobIdPrefixe = shift;

	my $trimJobIdVarNameSample = undef;
	my $libraryType = LoadConfig::getParam($rH_cfg, 'default', 'libraryType');
	for my $rH_laneInfo (@$rAoH_sampleLanes) {
		#print "mkdir -p metrics/$sampleName/output_jobs reads/$sampleName/output_jobs\n";
		##get raw read count
# 		my $inputFile = LoadConfig::getParam($rH_cfg, 'default', 'rawReadDir') .'/' .$sampleName .'/run' .$rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'} .'/' .$rH_laneInfo->{'read1File'};
# 		my $outputFile= 'metrics/' .$sampleName .'.' .$rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'} . '.readstats.raw.csv' ;
# 		my $command = Metrics::readStats($rH_cfg,$inputFile,$outputFile,'fastq',$libraryType);
# 		my $rawReadStatJobID = undef;
# 		if(defined($command) && length($command) > 0) {
# 			$rawReadStatJobID = SubmitToCluster::printSubmitCmd($rH_cfg, "metrics", 'raw', 'RAWREADSTAT' .$rH_jobIdPrefixe ->{$sampleName .'.' .$rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'}} , undef, $sampleName, $command);
# 			$rawReadStatJobID = '$'.$rawReadStatJobID;
# 		}
# 		
		## trimming - TO DO should be modified to the new rawread location (cf. David modif) and portability
		my $minQuality  = $rH_cfg->{'trim.minQuality'};
		my $minLength   = $rH_cfg->{'trim.minLength'};
		my $laneDir = 'reads/' .$sampleName . "/run" . $rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'} . "/";
		print "mkdir -p $laneDir\n";
		my $outputFastqPair1Name;
		if ( $rH_laneInfo->{'runType'} eq "SINGLE_END" ) {
			$outputFastqPair1Name = $laneDir . $sampleName.'.t'.$minQuality.'l'.$minLength.'.single.fastq.gz';
		}
		elsif ( $rH_laneInfo->{'runType'} eq "PAIRED_END" ) {
			$outputFastqPair1Name = $laneDir . $sampleName.'.t'.$minQuality.'l'.$minLength.'.pair1.fastq.gz';
		}
		else {
			die "Unknown runType: " . $rH_laneInfo->{' runType '} . "\n";
		}
		my $rH_trimDetails = Trimmomatic::trim($rH_cfg, $sampleName, $rH_laneInfo, $laneDir);
		my $trimJobIdVarNameLane=undef;
		if(length($rH_trimDetails->{'command'}) > 0) {
			$trimJobIdVarNameLane = SubmitToCluster::printSubmitCmd($rH_cfg, "trim", $rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'}, 'TRIM' .$rH_jobIdPrefixe ->{$sampleName.'.' .$rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'}} , undef, $sampleName, $rH_trimDetails->{'command'});
			$trimJobIdVarNameLane = '$' .$trimJobIdVarNameLane ;
			$trimJobIdVarNameSample .= $trimJobIdVarNameLane .LoadConfig::getParam($rH_cfg, 'default', 'clusterDependencySep');
		}
		###mbourgey 2013/08/01 - deprecated since louis parse directly the output of trimmomatic in the Trimmomatic.pm
		###new output stats file= 'read/' .$sampleName .'/' .$rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'} .'/'.$sampleName .'.trim.stats.csv'
# 			my $trinityOut = $laneDir .'/' . $sampleName . '.trim.out';
# 			##get trimmed read count
# 			my $outputFile= 'metrics/' .$sampleName .'.' .$rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'} . '.readstats.triming.tsv' ;
# 			my $command = Metrics::readStats($rH_cfg,$trinityOut,$outputFile,$sampleName.'.' .$rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'},'trim');
# 			my $filteredReadStatJobID ;
# 			if(defined($command) && length($command) > 0) {
# 				$filteredReadStatJobID = SubmitToCluster::printSubmitCmd($rH_cfg, "metrics", 'filtered', 'FILTERREADSTAT' .$rH_jobIdPrefixe ->{$sampleName.'.' .$rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'}} ,$trimJobIdVarNameLane, $sampleName, $command);
# 				$filteredReadStatJobID = '$'.$filteredReadStatJobID;
# 				$trimJobIdVarNameSample .= $filteredReadStatJobID .LoadConfig::getParam($rH_cfg, 'default', 'clusterDependencySep');
# 			}
	}
	if(defined($trimJobIdVarNameSample) && length($trimJobIdVarNameSample) > 0) {
		$trimJobIdVarNameSample = substr $trimJobIdVarNameSample, 0, -1 ;
	}
	return $trimJobIdVarNameSample;	
}


sub trimMetrics {
	my $depends = shift;
	my $rH_cfg = shift;
	my $rHoAoH_sampleInfo = shift;
	my $rHoAoA_designGroup  = shift;
	my $rAoH_seqDictionary = shift;
	my $rH_jobIdPrefixe = shift;

	my $libraryType = LoadConfig::getParam($rH_cfg, 'default', 'libraryType');
	my $trimmingDependency = undef;
	if($depends > 0) {
		$trimmingDependency = join(LoadConfig::getParam($rH_cfg, 'default', 'clusterDependencySep'),values(%{$globalDep{'trimming'}}));
	}
	print "mkdir -p metrics/\n";
	my $folder = 'reads';
	my $pattern = 'trim.stats.csv';
	my $ouputFile = 'metrics/trimming.stats';
	my $command;
	$command = Metrics::mergeTrimmomaticStats($rH_cfg,  $libraryType, $pattern, $folder, $ouputFile);
	my $metricsJobId = undef;
	if(defined($command) && length($command) > 0) {
		my $trimMetricsJobId = SubmitToCluster::printSubmitCmd($rH_cfg, "trimMetrics", undef, 'TRIMMETRICS', $trimmingDependency, undef, $command);
		$metricsJobId = '$' .$trimMetricsJobId;
	}
	return $metricsJobId;
}

sub aligning {
	my $depends = shift;
	my $rH_cfg = shift;
	my $sampleName = shift;
	my $rAoH_sampleLanes  = shift;
	my $rAoH_seqDictionary = shift;
	my $rH_jobIdPrefixe = shift;

	my $alignJobIdVarNameSample = undef;
	my $jobDependency = undef;
	if($depends > 0) {
	$jobDependency = $globalDep{'trimming'}{$sampleName};
	}
	
	for my $rH_laneInfo (@$rAoH_sampleLanes) {
		my $alignJobIdVarNameLane;
		my $pair1;
		my $pair2;
		my $single;
		my $command;
		#align lanes
		my $outputDirPath = 'alignment/' .$sampleName . "/run" . $rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'} . "/";
		print "mkdir -p $outputDirPath \n" ;
		if ( $rH_laneInfo->{'runType'} eq "SINGLE_END" ) {
			$single =  'reads/' .$sampleName . "/run" . $rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'} . "/" . $sampleName .'.t' .LoadConfig::getParam($rH_cfg,'trim','minQuality') .'l' .LoadConfig::getParam($rH_cfg,'trim','minLength') .'.single.fastq.gz';
			$command = TophatBowtie::align($rH_cfg, $sampleName, $rH_laneInfo, $single, ' ' );
		}
		elsif($rH_laneInfo->{'runType'} eq "PAIRED_END") {
			$pair1 =  'reads/' .$sampleName . "/run" . $rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'} . "/" . $sampleName .'.t' .LoadConfig::getParam($rH_cfg,'trim','minQuality') .'l' .LoadConfig::getParam($rH_cfg,'trim','minLength') .'.pair1.fastq.gz';
			$pair2 =  'reads/' .$sampleName . "/run" . $rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'} . "/" . $sampleName .'.t' .LoadConfig::getParam($rH_cfg,'trim','minQuality') .'l' .LoadConfig::getParam($rH_cfg,'trim','minLength') .'.pair2.fastq.gz';
			$command = TophatBowtie::align($rH_cfg, $sampleName, $rH_laneInfo, $pair1, $pair2);
		}
		if(defined $command && length($command) > 0){
			$alignJobIdVarNameLane = SubmitToCluster::printSubmitCmd($rH_cfg, "align", $rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'}, 'ALIGN' .$rH_jobIdPrefixe ->{$sampleName.'.' .$rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'}} .'ALIGN', $jobDependency, $sampleName, $command);
			$alignJobIdVarNameSample .= '$'. $alignJobIdVarNameLane .LoadConfig::getParam($rH_cfg, 'default', 'clusterDependencySep'); 
		} 
		##generate aligment stats
# 		my $inputFile = 'alignment/' . $sampleName . "/run" . $rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'} . "/" . 'accepted_hits.bam';
# 		my $outputFile= 'metrics/' .$sampleName .'.' .$rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'} . '.readstats.aligned.csv' ;
# 		$command = Metrics::readStats($rH_cfg,$inputFile,$outputFile,'bam');
# 		my $alignedReadStatJobID = undef;
# 		if(defined($command) && length($command) > 0) {
# 			$alignedReadStatJobID = SubmitToCluster::printSubmitCmd($rH_cfg, "metrics", 'aligned', 'ALIGNEDREADSTAT' .$rH_jobIdPrefixe ->{$sampleName.'.' .$rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'}} ,$alignJobIdVarNameLane, $sampleName, $command);
# 			$alignedReadStatJobID = '$'.$alignedReadStatJobID;
# 		}
# 		##merge read stats
# 		my $rawFile = 'metrics/' .$sampleName .'.' .$rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'} . '.readstats.raw.csv' ;
# 		my $filterFile = 'metrics/' .$sampleName .'.' .$rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'} . '.readstats.filtered.csv' ;
# 		my $alignFile = $outputFile ;
# 		my $sampleNameFull = $sampleName .'.' .$rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'};
# 		$outputFile= 'metrics/' .$sampleName .'.' .$rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'} . '.readstats.csv' ;
# 		$command = Metrics::mergeIndvidualReadStats($rH_cfg, $sampleNameFull, $rawFile, $filterFile, $alignFile, $outputFile);
# 		my $mergeReadStatJobID = undef;
# 		if(defined($command) && length($command) > 0) {
# 			$mergeReadStatJobID = SubmitToCluster::printSubmitCmd($rH_cfg, "metrics", 'merged', 'MERGEREADSTAT' .$rH_jobIdPrefixe ->{$sampleName.'.' .$rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'}} ,$alignedReadStatJobID, $sampleName, $command);
# 			$mergeReadStatJobID = '$'.$mergeReadStatJobID;
# 			$alignJobIdVarNameSample .= $mergeReadStatJobID .LoadConfig::getParam($rH_cfg, 'default', 'clusterDependencySep');
# 		}
		
	}
	$alignJobIdVarNameSample = substr $alignJobIdVarNameSample, 0, -1 ;
	return $alignJobIdVarNameSample;
}

sub merging {
	my $depends = shift;
	my $rH_cfg = shift;
	my $sampleName = shift;
	my $rAoH_sampleLanes  = shift;
	my $rAoH_seqDictionary = shift;
	my $rH_jobIdPrefixe = shift;

	my $jobDependency = undef;
	if($depends > 0) {
		$jobDependency = $globalDep{'aligning'}{$sampleName};
	}
	
	##Merging
	my $inputBAM ; 
	my $outputBAM = 'alignment/' .$sampleName .'/' .$sampleName .'.merged.bam' ;
	my @alignFiles;
	my $merge = 0; # JT : Flag to see if we skip merging step or not.

	for my $rH_laneInfo (@$rAoH_sampleLanes) {
		my $laneDir = "alignment/" . $sampleName . "/run" . $rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'} . "/";
		$inputBAM = $laneDir . 'accepted_hits.bam';
		$merge++; # JT: increment flag 
		push(@alignFiles, $inputBAM) ;
	}
	my $command;
	my $mergeJobId = undef;
	if($merge > 1){ # JT: If flag is higher than 1 (so more than one lane / sample), perform merge.
		$command  = Picard::mergeFiles($rH_cfg, $sampleName, \@alignFiles, $outputBAM);
		if(defined($command) && length($command) > 0) {
			$mergeJobId = SubmitToCluster::printSubmitCmd($rH_cfg, "mergeFiles", undef, 'MERGELANES' .$rH_jobIdPrefixe ->{$sampleName} , $jobDependency, $sampleName, $command);
			$mergeJobId = '$'.$mergeJobId;
		}
	}else{
		$mergeJobId = $jobDependency; # JT : update job dependency as well.
	}
	
	## reorder
	if($merge > 1){ # JT : if merge, update input file accordingly.
		$inputBAM = $outputBAM; #Update file name if it has been merged.
	}
	$outputBAM = 'alignment/' .$sampleName .'/' .$sampleName .'.merged.karyotypic.bam';
	$command = Picard::reorderSam($rH_cfg, $sampleName, $inputBAM, $outputKarBAM);
	my $reorderJobId = undef;
	if(defined($command) && length($command) > 0) {
		$reorderJobId = SubmitToCluster::printSubmitCmd($rH_cfg, "reorderSam", undef, 'REORDER' .$rH_jobIdPrefixe ->{$sampleName} .'REORDER', $mergeJobId, $sampleName, $command);
		$reorderJobId = '$'.$reorderJobId;
	}
	## mark duplicates
	$inputBAM = $outputBAM;
	$outputBAM = 'alignment/' .$sampleName .'/' .$sampleName .'.merged.mdup.bam';
	my $duplicatesMetricsFile = 'alignment/' .$sampleName .'/' .$sampleName .'.merged.mdup.metrics';
	$command = Picard::markDup($rH_cfg, $sampleName, $inputBAM, $outputBAM,$duplicatesMetricsFile );
	my $markDupJobId = undef;
	if(defined($command) && length($command) > 0) {
		$markDupJobId = SubmitToCluster::printSubmitCmd($rH_cfg, "markDup", undef, 'MARKDUP' .$rH_jobIdPrefixe ->{$sampleName} , $reorderJobId, $sampleName, $command);
		$markDupJobId = '$'.$markDupJobId
	}
	return $markDupJobId;
}



sub alignMetrics {
	my $depends = shift;
	my $rH_cfg = shift;
	my $rHoAoH_sampleInfo = shift;
	my $rHoAoA_designGroup  = shift;
	my $rAoH_seqDictionary = shift;
	my $rH_jobIdPrefixe = shift;

	my $mergingDependency = undef;
	if($depends > 0) {
		$mergingDependency = join(LoadConfig::getParam($rH_cfg, 'default', 'clusterDependencySep'),values(%{$globalDep{'merging'}}));
	}
	## RNAseQC metrics
	mkdir  $workDir .'/alignment' ;
	my $sampleList = $workDir .'/alignment/rnaseqc.samples.txt';
	open(RNASAMPLE, ">$sampleList") or  die ("Unable to open $sampleList for writing") ;
	print RNASAMPLE "Sample\tBamFile\tNote\n";
	my $projectName = LoadConfig::getParam($rH_cfg, 'metricsRNA', 'projectName');
	for my $sampleName (keys %{$rHoAoH_sampleInfo}) {
		print RNASAMPLE "$sampleName\talignment/$sampleName/$sampleName.merged.mdup.bam\t$projectName\n";
	}
	print "mkdir -p metrics/rnaseqRep/\n";
	$sampleList = 'alignment/rnaseqc.samples.txt';
	my $outputFolder = 'metrics/rnaseqRep';
	my $command = Metrics::rnaQc($rH_cfg, $sampleList, $outputFolder);
	my $rnaqcJobId = undef;
	my $metricsJobId = undef;
	if(defined($command) && length($command) > 0) {
		$rnaqcJobId = SubmitToCluster::printSubmitCmd($rH_cfg, "rnaQc", undef, 'METRICSRNA', $mergingDependency, undef, $command);
		$metricsJobId .= '$' .$rnaqcJobId ;
	}

	#$metricsJobId = substr $metricsJobId, 0, -1 ;
	return $metricsJobId;
}


# sub mutation{
# 	my $depends = shift;
# 	my $rH_cfg = shift;
# 	my $sampleName = shift;
# 	my $rAoH_sampleLanes  = shift;
# 	my $rAoH_seqDictionary = shift;
# ALIGNEDREADSTAT_JOB_ID
# 	my $jobDependency = undef;
# 	if($depends > 0) {
# 		$jobDependency = $globalDep{'merging'}{$sampleName};
# 	}
# 
# 	my $inputBAM = "alignment/" . $sampleName . "/" . $sampleName . ".merged.mdup.bam" ; 
# 
# 
# }

sub wiggle {
	my $depends = shift;
	my $rH_cfg = shift;
	my $sampleName = shift;
	my $rAoH_sampleLanes  = shift;
	my $rAoH_seqDictionary = shift;
	my $rH_jobIdPrefixe = shift;

	my $jobDependency = undef;
	if($depends > 0) {
		$jobDependency = $globalDep{'merging'}{$sampleName};
	}

	my $inputBAM = "alignment/" . $sampleName . "/" . $sampleName . ".merged.mdup.bam" ; 
	#testing for strand-specificity

	my $strandSPecificityInfo = LoadConfig::getParam($rH_cfg, 'align', 'strandInfo');
	my @strandJobId ;
	my @outputBAM;
	my @outputBedGraph;
	my @outputWiggle;
	my @prefixJobName;
	print "mkdir -p tracks/$sampleName/ tracks/bigWig/\n";
	if($strandSPecificityInfo ne "fr-unstranded") {
	## strand specific 
		@outputBAM = ('alignment/' . $sampleName . '/' . $sampleName . '.merged.mdup.forward.bam' ,  'alignment/' . $sampleName . '/' . $sampleName . '.merged.mdup.reverse.bam');
		my $rA_command = Wiggle::strandBam($rH_cfg, $sampleName, $inputBAM, \@outputBAM);
		if(defined($rA_command) && @{$rA_command} > 1) {
			my $strandJobIdF = SubmitToCluster::printSubmitCmd($rH_cfg, "wiggle", 'FORWARD1', 'FSTRANDSPEC' .$rH_jobIdPrefixe ->{$sampleName} , $jobDependency, $sampleName, $rA_command->[0]);
			push(@strandJobId, '$'.$strandJobIdF );
			my $strandJobIdR = SubmitToCluster::printSubmitCmd($rH_cfg, "wiggle", 'REVERSE1', 'RSTRANDSPEC' .$rH_jobIdPrefixe ->{$sampleName} , $jobDependency, $sampleName, $rA_command->[1]);
			push(@strandJobId, '$'.$strandJobIdR );
		}
		@outputBedGraph = ('tracks/' . $sampleName . '/' . $sampleName . '.forward.bedGraph' ,  'tracks/' . $sampleName . '/' . $sampleName . '.reverse.bedGraph');
		@outputWiggle = ('tracks/bigWig/' . $sampleName . '.forward.bw' ,  'tracks/bigWig/' . $sampleName . '.reverse.bw');
		@prefixJobName = ( 'FORWARD2', 'REVERSE2');
	}
	else {
		push(@outputBAM,$inputBAM);
		push(@strandJobId, $jobDependency);
		push(@outputBedGraph,'tracks/' . $sampleName . '/' . $sampleName . '.bedGraph');
		push(@outputWiggle,'tracks/bigWig/' . $sampleName . '.bw' );
		push(@prefixJobName , undef ) ;
	}
	my $wiggleJobId ;
	for(my $i = 0; $i <@outputBAM; $i++) {
		my $command = Wiggle::graph($rH_cfg, $sampleName, $inputBAM, $outputBedGraph[$i], $outputWiggle[$i]);
		if(defined($command) && length($command) > 0) {
			my $tmpJobId = SubmitToCluster::printSubmitCmd($rH_cfg, "wiggle", $prefixJobName[$i], 'WIGGLE' .$rH_jobIdPrefixe ->{$sampleName} , $strandJobId[$i], $sampleName, $command);
			$wiggleJobId .= '$'.$tmpJobId .LoadConfig::getParam($rH_cfg, 'default', 'clusterDependencySep');
		}
	} 
	$wiggleJobId = substr $wiggleJobId, 0, -1 ;
	return $wiggleJobId;	
}


sub rawCounts {
	my $depends = shift;
	my $rH_cfg = shift;
	my $sampleName = shift;
	my $rAoH_sampleLanes  = shift;
	my $rAoH_seqDictionary = shift;
	my $rH_jobIdPrefixe = shift;

	my $jobDependency = undef;
	if($depends > 0) {
		$jobDependency = $globalDep{'merging'}{$sampleName};
	}
	print "mkdir -p raw_counts\n";
	my $inputBAM = 'alignment/' . $sampleName . '/' . $sampleName . '.merged.mdup.bam' ;
        my $sortedBAM = 'alignment/' . $sampleName . '/' . $sampleName . '.queryNameSorted.bam' ;
	my $inputGtf = LoadConfig::getParam($rH_cfg, 'htseq', 'referenceGtf');
	my $outputCount = 'raw_counts/' . $sampleName . '.readcounts.csv';
	my $sortOrder = 'queryname';
	my $strandInfo;
	my $strandSPecificityInfo = LoadConfig::getParam($rH_cfg, 'align', 'strandInfo');
	if($strandSPecificityInfo ne "fr-unstranded") {
		 $strandInfo= 'reverse';
	}
	else {
		$strandInfo= 'no';
	}
	## query sort the bam
	my $sortJobId;
	my $command = Picard::sortSam($rH_cfg, $sampleName, $inputBAM, $sortedBAM, $sortOrder);
	if(defined($command) && length($command) > 0) {
		$sortJobId=SubmitToCluster::printSubmitCmd($rH_cfg, "sortSam", undef, 'QNSORT' .$rH_jobIdPrefixe ->{$sampleName} , $jobDependency, $sampleName, $command);
		$sortJobId='$'.$sortJobId
	}
	## count reads
        my $countJobId;
	$command = HtseqCount::readCountPortable($rH_cfg, $sortedBAM, $inputGtf, $outputCount, $strandInfo); 
	if(defined($command) && length($command) > 0) {
		$countJobId=SubmitToCluster::printSubmitCmd($rH_cfg, "htseq", undef, 'RAWCOUNT' .$rH_jobIdPrefixe ->{$sampleName} , $sortJobId, $sampleName, $command);
		$countJobId='$'.$countJobId
	}
	return $countJobId;
}


sub rawCountsMetrics {
  my $depends = shift;
  my $rH_cfg = shift;
  my $rHoAoH_sampleInfo = shift;
  my $rHoAoA_designGroup  = shift;
  my $rAoH_seqDictionary = shift;
  my $rH_jobIdPrefixe = shift;

	my $metricsJobId  = undef;
	my $countDependency = undef;
	my $wiggleDependency = undef;
	if($depends > 0) {
		$countDependency = join(LoadConfig::getParam($rH_cfg, 'default', 'clusterDependencySep'),values(%{$globalDep{'rawCounts'}}));
		$wiggleDependency = join(LoadConfig::getParam($rH_cfg, 'default', 'clusterDependencySep'),values(%{$globalDep{'wiggle'}}));
	}
	#create rawcount matrix
	print "mkdir -p DGE\n";
	my $readCountDir = 'raw_counts' ;
	my $readcountExtension = '.readcounts.csv';
	my $outputDir = 'DGE';
	my $outputMatrix = 'rawCountMatrix.csv';
	my $command = HtseqCount::refGtf2matrix($rH_cfg, LoadConfig::getParam($rH_cfg, 'htseq', 'referenceGtf'), $readCountDir, $readcountExtension, $outputDir, $outputMatrix);
	my $matrixJobId = undef;
	if(defined($command) && length($command) > 0) {
		$matrixJobId = SubmitToCluster::printSubmitCmd($rH_cfg, "metrics", 'matrix', 'MATRIX', $countDependency, undef, $command);
		$matrixJobId = '$' .$matrixJobId;
		$metricsJobId .= $matrixJobId .LoadConfig::getParam($rH_cfg, 'default', 'clusterDependencySep');
	}
	
	### to do outside of the wiggle function on time only
	print "mkdir -p tracks/bigWig\n";
	my $wigFolder = 'tracks/bigWig/' ;
	my $wigArchive = 'tracks.zip' ;
	$command = Wiggle::zipWig($rH_cfg, $wigFolder, $wigArchive);
	my $wigZipJobId ;
	if(defined($command) && length($command) > 0) {
	    my $tmpJobId = SubmitToCluster::printSubmitCmd($rH_cfg, "metrics", undef, 'WIGZIP' , $wiggleDependency, undef, $command);
	    $metricsJobId .= '$' .$tmpJobId .LoadConfig::getParam($rH_cfg, 'default', 'clusterDependencySep');
	}
	##RPKM and Saturation
	print "mkdir -p metrics/saturation\n";;
	my $countFile   = 'DGE/rawCountMatrix.csv';
	my $geneSizeFile     = LoadConfig::getParam($rH_cfg, 'saturation', 'geneSizeFile');
	my $rpkmDir = 'raw_counts';
	my $saturationDir = 'metrics/saturation';
	
	$command =  Metrics::saturation($rH_cfg, $countFile, $geneSizeFile, $rpkmDir, $saturationDir);
	my $saturationJobId = undef;
	if(defined($command) && length($command) > 0) {
		$saturationJobId = SubmitToCluster::printSubmitCmd($rH_cfg, "saturation", undef, 'RPKM', $matrixJobId, undef, $command);
		$metricsJobId .= '$' .$saturationJobId .LoadConfig::getParam($rH_cfg, 'default', 'clusterDependencySep');
	}
	$metricsJobId = substr $metricsJobId, 0, -1 ;
	return $metricsJobId;
}

sub fpkm {
	my $depends = shift;
	my $rH_cfg = shift;
	my $sampleName = shift;
	my $rAoH_sampleLanes  = shift;
	my $rAoH_seqDictionary = shift;
	my $rH_jobIdPrefixe = shift;

	my $jobDependency = undef;
	if($depends > 0) {
		$jobDependency = $globalDep{'merging'}{$sampleName};
	}
	print "mkdir -p fpkm/known/$sampleName fpkm/denovo/$sampleName\n";
	my $inputBAM = 'alignment/' . $sampleName . '/' . $sampleName . '.merged.mdup.bam' ;
	my $outputKnown = 'fpkm/known/' . $sampleName;
	my $outputDeNovo = 'fpkm/denovo/' . $sampleName;
	my $gtfOption = '-G ' .LoadConfig::getParam($rH_cfg, 'fpkm','referenceGtf');
	
	
	## known FPKM
	my $fpkmJobId = undef;
	my $command = Cufflinks::fpkm($rH_cfg, $inputBAM, $outputKnown, $gtfOption);
	if(defined($command) && length($command) > 0) {
		my $fpkmKnownJobId=SubmitToCluster::printSubmitCmd($rH_cfg, "fpkm", "KNOWN", 'FPKMK' .$rH_jobIdPrefixe ->{$sampleName} .'FPKM', $jobDependency, $sampleName, $command);
		$fpkmJobId = '$' .$fpkmKnownJobId .LoadConfig::getParam($rH_cfg, 'default', 'clusterDependencySep');
	}
	## denovo FPKM
	$command = Cufflinks::fpkm($rH_cfg, $inputBAM, $outputDeNovo, undef);
	if(defined($command) && length($command) > 0) {
		my $fpkmDeNovoJobId=SubmitToCluster::printSubmitCmd($rH_cfg, "fpkm", "DENOVO", 'FPKMD' .$rH_jobIdPrefixe ->{$sampleName} , $jobDependency, $sampleName, $command);
		$fpkmJobId .= '$' .$fpkmDeNovoJobId  .LoadConfig::getParam($rH_cfg, 'default', 'clusterDependencySep');
	}
	if (defined($fpkmJobId)) {
		$fpkmJobId = substr $fpkmJobId, 0, -1 ;
	}
	return $fpkmJobId;
}


sub cuffdiff {
	my $depends = shift;
	my $rH_cfg = shift;
	my $rHoAoH_sampleInfo = shift;
	my $rHoAoA_designGroup  = shift;
	my $rAoH_seqDictionary = shift;
	my $rH_jobIdPrefixe = shift;
	
	my $jobDependency = undef;
	my $cuffJobID ;
	if($depends > 0 and values(%{$globalDep{'fpkm'}}) > 0) {
		$jobDependency = join(LoadConfig::getParam($rH_cfg, 'default', 'clusterDependencySep'),values(%{$globalDep{'fpkm'}}));
	}
	print "mkdir -p cuffdiff/known cuffdiff/denovo\n";
	## create the list of deNovo gtf to merge
	mkdir $workDir .'/fpkm';
	mkdir $workDir .'/fpkm/denovo/';
	my $mergeListFile = $workDir .'/fpkm/denovo/gtfMerge.list';
	my $compareList = " ";
	open(MERGEF, ">$mergeListFile") or  die ("Unable to open $mergeListFile for wrtting") ;
	##iterate over sample
	for my $sampleName (keys %{$rHoAoH_sampleInfo}) {
	    my $gtfFile = 'fpkm/denovo/' .$sampleName .'/transcripts.gtf' ;
	    $compareList .= 'fpkm/denovo/' .$sampleName .'/transcripts.gtf ' ;
 	    print MERGEF $gtfFile;
	    print MERGEF "\n";
	}
	close($mergeListFile);
	##merge denovo transcript in one  gtf file
 	my $outputPathDeNovo = 'fpkm/denovo/allSample' ;
 	my $command = Cufflinks::cuffcompare($rH_cfg, $compareList, $outputPathDeNovo, $mergeListFile);
 	my $cuffmergeJobId ;
 	if(defined($command) && length($command) > 0) {
 	    $cuffmergeJobId = SubmitToCluster::printSubmitCmd($rH_cfg, "cuffcompare", "MERGE", 'GTFCOMPARE', $jobDependency, undef, $command);
 	    $cuffmergeJobId = '$' .$cuffmergeJobId 
 	}
 	my $gtfDnMerged = 'fpkm/denovo/merged.gtf';
 	my $gtfDnFormatMerged = 'fpkm/denovo/formated.merged.gtf';
#	$command = Cufflinks::mergeGtfFormat($rH_cfg, $gtfDnMerged, $gtfDnFormatMerged);
# 	my $formatJobId;
# 	if(defined($command) && length($command) > 0) {
# 	    $formatJobId= SubmitToCluster::printSubmitCmd($rH_cfg, "default", "FORMAT", 'GTFFORMAT', $cuffmergeJobId, undef, $command);
# 	    $formatJobId= '$' .$formatJobId;
#	    $cuffJobID = $formatJobId .LoadConfig::getParam($rH_cfg, 'default', 'clusterDependencySep');
# 	}
	##iterate over design
	my $cuffddiffJobId;
	for my $design (keys %{$rHoAoA_designGroup}) {
		mkdir $workDir;
		mkdir $workDir .'/cuffdiff';
		my $numberGroups = @{$rHoAoA_designGroup->{$design}} ;
		##iterate over group
		my @groupInuptFiles;
		for (my $i = 0;   $i < $numberGroups; $i++) {
			##iterate over samples in the design
			my $numberSample =  @{$rHoAoA_designGroup->{$design}->[$i]};
			my $bamfile = ' ';
			for (my $j = 0;   $j < $numberSample; $j++) {
				$bamfile .= 'alignment/' .$rHoAoA_designGroup->{$design}->[$i]->[$j] . '/' .$rHoAoA_designGroup->{$design}->[$i]->[$j] . '.merged.mdup.bam' .',' ;
			}
			$bamfile = substr $bamfile, 0, -1 ;
			push(@groupInuptFiles,$bamfile);
		}

		my $outputPathKnown = 'cuffdiff/known/' .$design;
		##cuffdiff known
		$command = Cufflinks::cuffdiff($rH_cfg,\@groupInuptFiles,$outputPathKnown,LoadConfig::getParam($rH_cfg, 'cuffdiff','referenceGtf'));
		if(defined($command) && length($command) > 0) {
			my $diffJobId = SubmitToCluster::printSubmitCmd($rH_cfg, "cuffdiff", "KNOWN",  'CUFFDIFFK' .$rH_jobIdPrefixe ->{$design} , $jobDependency, $design, $command);
			$cuffddiffJobId .= '$' .$diffJobId .LoadConfig::getParam($rH_cfg, 'default', 'clusterDependencySep');
		}
	}
	if(defined($cuffddiffJobId) && length($cuffddiffJobId) > 0) {
		$cuffddiffJobId = substr $cuffddiffJobId, 0, -1 ;
	}
	$command = Cufflinks::mergeCuffdiffRes($rH_cfg,$designFilePath,'cuffdiff','fpkm');
	my $mergeCuffdiffResJobID;
	if(defined($command) && length($command) > 0) {
		$mergeCuffdiffResJobID = SubmitToCluster::printSubmitCmd($rH_cfg, "default", "MERGE_RES", 'CUFF_MERGE_RES', $cuffddiffJobId, undef, $command);
		$mergeCuffdiffResJobID = '$' .$mergeCuffdiffResJobID;
		$cuffJobID .= $mergeCuffdiffResJobID .LoadConfig::getParam($rH_cfg, 'default', 'clusterDependencySep');
	}
	### mbourgey 29/07/2013 - filtering now included in the R script that merge cuffdiff res with fpkm - remove old filtering 
	if(defined($cuffJobID) && length($cuffJobID) > 0) {
	    $cuffJobID = substr $cuffJobID, 0, -1 ;
	}
	return $cuffJobID;
}


sub dge {
	my $depends = shift;
	my $rH_cfg = shift;
	my $rHoAoH_sampleInfo = shift;
	my $rHoAoA_designGroup  = shift;
	my $rAoH_seqDictionary = shift;
	my $rH_jobIdPrefixe = shift;

	my $jobDependency = undef;
	if($depends > 0) {
		$jobDependency = $globalDep{'rawCountsMetrics'}{'rawCountsMetrics'};
	}
	
	my $countMatrix = 'DGE/rawCountMatrix.csv';
	my $outputDir = 'DGE';
	
	## edgeR
	my $command = DiffExpression::edgerPortable($rH_cfg, $designFilePath, $countMatrix, $outputDir);
	my $edgerJobId = undef;
	if(defined($command) && length($command) > 0) {
		$edgerJobId = SubmitToCluster::printSubmitCmd($rH_cfg, "diffExpress", 'edger', 'EDGER', $jobDependency, undef, $command);
		$edgerJobId = '$' .$edgerJobId;
	}
	
	## DESeq
	$command = DiffExpression::deseq($rH_cfg, $designFilePath, $countMatrix, $outputDir);
	my $deseqJobId = undef;
	if(defined($command) && length($command) > 0) {
		$deseqJobId = SubmitToCluster::printSubmitCmd($rH_cfg, "diffExpress", 'deseq', 'DESEQ', $edgerJobId, undef,$command);
		$deseqJobId = '$' .$deseqJobId;
	}
	
	return $deseqJobId;
} 

sub goseq {
	my $depends = shift;
	my $rH_cfg = shift;
	my $rHoAoH_sampleInfo = shift;
	my $rHoAoA_designGroup  = shift;
	my $rAoH_seqDictionary = shift;
	my $rH_jobIdPrefixe = shift;

	
	my $cuffdiffDependency = undef;
	my $dgeDependency = undef;
	if($depends > 0) {
		$cuffdiffDependency = $globalDep{'cuffdiff'}{'cuffdiff'};
		$dgeDependency = $globalDep{'dge'}{'dge'};
	}

	my $columnsCuff = LoadConfig::getParam($rH_cfg, 'diffExpress', 'cuffRescolumns');
	my $columnsDge = LoadConfig::getParam($rH_cfg, 'diffExpress', 'dgeRescolumns');
	my $command;
	my $goseqJobId;
	for my $design (keys %{$rHoAoA_designGroup}) {
		## goseq for cuffdiff known results
		my $resultFileCuff = 'cuffdiff/known/' .$design .'/isoform_exp.diff' ;
		my $outputFileCuff = 'cuffdiff/known/' .$design .'/gene_ontology_results.csv';
		$command = DiffExpression::goseq($rH_cfg, $resultFileCuff, $outputFileCuff, $columnsCuff);
		my $goCuffJobId = undef;
		if(defined($command) && length($command) > 0) {
			$goCuffJobId = SubmitToCluster::printSubmitCmd($rH_cfg, "goseq","CUFFLINKS", 'GOCUFFDIFF' .$rH_jobIdPrefixe ->{$design} , $cuffdiffDependency, $design, $command);
			$goseqJobId .= '$' .$goCuffJobId .LoadConfig::getParam($rH_cfg, 'default', 'clusterDependencySep');
		}
		## goseq for dge results
		my $resultFileDge = 'DGE/' .$design .'/dge_results.csv' ;
		my $outputFileDge = 'DGE/' .$design .'/gene_ontology_results.csv';
		$command = DiffExpression::goseq($rH_cfg, $resultFileDge, $outputFileDge, $columnsDge);
		my $goDgeJobId = undef;
		if(defined($command) && length($command) > 0) {
			$goDgeJobId = SubmitToCluster::printSubmitCmd($rH_cfg, "goseq", "DESEQ", 'GODGE' .$rH_jobIdPrefixe ->{$design} , $dgeDependency, $design, $command);
			$goseqJobId .= '$' .$goDgeJobId .LoadConfig::getParam($rH_cfg, 'default', 'clusterDependencySep');
		}
	}
	if(defined($goseqJobId) && length($goseqJobId) > 0) {
		$goseqJobId = substr $goseqJobId, 0, -1 ;
	}
	return $goseqJobId;
}

sub exploratory {
	my $depends = shift;
	my $rH_cfg = shift;
	my $rHoAoH_sampleInfo = shift;
	my $rHoAoA_designGroup  = shift;
	my $rAoH_seqDictionary = shift;
	my $rH_jobIdPrefixe = shift;

	my $jobDependency = undef;
	if($depends > 0 and values(%{$globalDep{'fpkm'}}) > 0) {
		$jobDependency   = join(LoadConfig::getParam($rH_cfg, 'default', 'clusterDependencySep'),values(%{$globalDep{'fpkm'}}));
		$jobDependency  .= LoadConfig::getParam($rH_cfg, 'default', 'clusterDependencySep') .$globalDep{'rawCountsMetrics'}{'rawCountsMetrics'} ;
	}

	print "mkdir -p exploratory\n";
	
	# Call gqSeqUtils::exploratoryAnalysis()
	# sub exploratoryRnaAnalysis{
	#         my $rH_cfg        = shift;
	#         my $readSetSheet  = shift;
	#         my $workDir       = shift;
	#         my $configFile    = shift;
	#
	my $command = GqSeqUtils::exploratoryRnaAnalysis($rH_cfg, $readSetSheet, $workDir, $configFile) ;
	my $exploratoryJobId = undef;
	if(defined($command) && length($command) > 0) {
		$exploratoryJobId = SubmitToCluster::printSubmitCmd($rH_cfg, "exploratory", 'exploratoryAnalysis', 'exploratoryAnalysis' , $jobDependency ,undef, $command);
		$exploratoryJobId = '$' .$exploratoryJobId ;
	}
	
	return $exploratoryJobId;
} 

sub deliverable {
	my $depends = shift;
	my $rH_cfg = shift;
	my $rHoAoH_sampleInfo = shift;
	my $rHoAoA_designGroup   = shift;
	my $rAoH_seqDictionary = shift;
	my $rH_jobIdPrefixe = shift;

	my $jobDependency ;
	if($depends > 0) {
		my $trimDependency = $globalDep{'trimMetrics'}{'trimMetrics'};
                if (defined($trimDependency) && length($trimDependency) > 0) {
                         $jobDependency .= $trimDependency .LoadConfig::getParam($rH_cfg, 'default', 'clusterDependencySep');
                }
		my $alignDependency = $globalDep{'alignMetrics'}{'alignMetrics'};
		if (defined($alignDependency) && length($alignDependency) > 0) {
                         $jobDependency .= $alignDependency .LoadConfig::getParam($rH_cfg, 'default', 'clusterDependencySep');
                }
		my $goDependency = $globalDep{'goseq'}{'goseq'};
		if (defined($goDependency) && length($goDependency) > 0) {
			 $jobDependency .= $goDependency .LoadConfig::getParam($rH_cfg, 'default', 'clusterDependencySep');
		}
	}

	if (defined($jobDependency) && length($jobDependency) > 0) {
                 $jobDependency = substr $jobDependency, 0, -1 ;
        }

	my $command = GqSeqUtils::clientReport($rH_cfg,  $configFile, $workDir, 'RNAseq') ;

	my $deliverableJobId = undef;
        if(defined($command) && length($command) > 0) {
		print "mkdir -p deliverable\n";
                $deliverableJobId = SubmitToCluster::printSubmitCmd($rH_cfg, "deliverable", 'REPORT', 'RNAREPORT', $jobDependency , undef, $command);
                $deliverableJobId = '$' .$deliverableJobId ;
        }

        return $deliverableJobId;
}
1;
