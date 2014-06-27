#!/usr/bin/env python

# Python Standard Modules

# MUGQIC Modules
from core.config import *
from core.job import *

def base_recalibrator(input, output):

    job = Job([input], [output], [['baseRecalibrator', 'moduleVersion.java'], ['baseRecalibrator', 'moduleVersion.gatk']])

    job.command = \
"""java -Djava.io.tmpdir={tmp_dir} {extra_java_flags} -Xmx{ram} -jar \$GATK_JAR \\
  -T BaseRecalibrator \\
  -nct {threads} \\
  -I {input} \\
  -R {reference_fasta} \\
  -knownSites {known_sites} \\
  -o {output}""".format(
        tmp_dir=config.param('baseRecalibrator', 'tmpDir'),
        extra_java_flags=config.param('baseRecalibrator', 'extraJavaFlags'),
        ram=config.param('baseRecalibrator', 'ram'),
        threads=config.param('baseRecalibrator', 'threads', type='int'),
        input=input,
        reference_fasta=config.param('baseRecalibrator', 'referenceFasta', type='filepath'),
        known_sites=config.param('baseRecalibrator', 'knownSites', type='filepath'),
        output=output
    )

    return job

def print_reads(input, output, base_quality_score_recalibration):

    job = Job([input], [output], [['printReads', 'moduleVersion.java'], ['printReads', 'moduleVersion.gatk']])

    job.command = \
"""java -Djava.io.tmpdir={tmp_dir} {extra_java_flags} -Xmx{ram} -jar \$GATK_JAR \\
  -T PrintReads \\
  -nct {threads} \\
  -I {input} \\
  -R {reference_fasta} \\
  -BQSR {base_quality_score_recalibration} \\
  -o {output}""".format(
        tmp_dir=config.param('printReads', 'tmpDir'),
        extra_java_flags=config.param('printReads', 'extraJavaFlags'),
        ram=config.param('printReads', 'ram'),
        threads=config.param('printReads', 'threads', type='int'),
        input=input,
        reference_fasta=config.param('printReads', 'referenceFasta', type='filepath'),
        base_quality_score_recalibration=base_quality_score_recalibration,
        output=output
    )

    return job


def realigner_target_creator(input, output, intervals=[], exclude_intervals=[]):

    job = Job([input], [output], [['realignerTargetCreator', 'moduleVersion.java'], ['realignerTargetCreator', 'moduleVersion.gatk']])

    job.command = \
"""java -Djava.io.tmpdir={tmp_dir} {extra_java_flags} -Xmx{ram} -jar \$GATK_JAR \\
  -T RealignerTargetCreator {extra_realigner_target_creator_flags} \\
  -R {reference_fasta} \\
  -I {input} \\
  -o {output}{intervals}{exclude_intervals}""".format(
        tmp_dir=config.param('realignerTargetCreator', 'tmpDir'),
        extra_java_flags=config.param('realignerTargetCreator', 'extraJavaFlags'),
        ram=config.param('realignerTargetCreator', 'ram'),
        extra_realigner_target_creator_flags=config.param('realignerTargetCreator', 'extraRealignerTargetCreatorFlags'),
        reference_fasta=config.param('realignerTargetCreator', 'referenceFasta', type='filepath'),
        input=input,
        output=output,
        intervals="".join(" \\\n  --intervals " + interval for interval in intervals),
        exclude_intervals="".join(" \\\n  --excludeIntervals " + exclude_interval for exclude_interval in exclude_intervals)
    )

    return job


def indel_realigner(input, output, target_intervals, intervals=[], exclude_intervals=[]):

    job = Job([input], [output], [['indelRealigner', 'moduleVersion.java'], ['indelRealigner', 'moduleVersion.gatk']])

    job.command = \
"""java -Djava.io.tmpdir={tmp_dir} {extra_java_flags} -Xmx{ram} -jar \$GATK_JAR \\
  -T IndelRealigner {extra_indel_realigner_flags} \\
  -R {reference_fasta} \\
  -I {input} \\
  --targetIntervals {target_intervals} \\
  -o {output}{intervals}{exclude_intervals}{extra_indel_realigner_flags} \\
  --maxReadsInMemory {max_reads_in_memory}""".format(
        tmp_dir=config.param('indelRealigner', 'tmpDir'),
        extra_java_flags=config.param('indelRealigner', 'extraJavaFlags'),
        ram=config.param('indelRealigner', 'ram'),
        extra_indel_realigner_flags=config.param('indelRealigner', 'extraIndelRealignerFlags'),
        reference_fasta=config.param('indelRealigner', 'referenceFasta', type='filepath'),
        input=input,
        target_intervals=target_intervals,
        output=output,
        intervals="".join(" \\\n  --intervals " + interval for interval in intervals),
        exclude_intervals="".join(" \\\n  --excludeIntervals " + exclude_interval for exclude_interval in exclude_intervals),
        max_reads_in_memory=config.param('indelRealigner', 'maxReadsInMemory')
    )

    return job

def depth_of_coverage(input, output_prefix, intervals=[]):

    job = Job([input], [output_prefix + ".sample_summary"], [['depthOfCoverage', 'moduleVersion.java'], ['depthOfCoverage', 'moduleVersion.gatk']])

    summary_coverage_thresholds = config.param('depthOfCoverage', 'percentThresholds', type='list').sort()

    job.command = \
"""java -Djava.io.tmpdir={tmp_dir} {extra_java_flags} -Xmx{ram} -jar \$GATK_JAR \\
  -T DepthOfCoverage --omitDepthOutputAtEachBase --logging_level ERROR \\
  -R {reference_fasta} \\
  -I {input} \\
  -o {output_prefix}{intervals}{summary_coverage_thresholds} \\
  --start 1 --stop {highest_summary_coverage_threshold} --nBins {highest_summary_coverage_threshold} -dt NONE""".format(
        tmp_dir=config.param('depthOfCoverage', 'tmpDir'),
        extra_java_flags=config.param('depthOfCoverage', 'extraJavaFlags'),
        ram=config.param('depthOfCoverage', 'ram'),
        reference_fasta=config.param('depthOfCoverage', 'referenceFasta', type='filepath'),
        input=input,
        output_prefix=output_prefix,
        intervals=" \\\n  --intervals ".join(intervals),
        summary_coverage_thresholds=" \\\n  --summaryCoverageThreshold ".join(summary_coverage_thresholds),
        highest_summary_coverage_threshold=summary_coverage_thresholds[-1]
    )

    return job
