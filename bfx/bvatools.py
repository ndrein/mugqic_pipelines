#!/usr/bin/env python

# Python Standard Modules

# MUGQIC Modules
from core.config import *
from core.job import *

def resolve_readset_coverage_bed(readset):
    coverage_target = config.param('bvatools_depth_of_coverage', 'coverage_targets', required=False)

    if coverage_target:
        if coverage_target == "auto":
            if readset.beds:
                return readset.beds[0]
            else:
                return None
        else:
            # Add filepath validation
            coverage_target = config.param('bvatools_depth_of_coverage', 'coverage_targets', type='filepath')
            return coverage_target
    else:
        return None

# If per RG != 0 is given there will be multiple outputs, so output is a prefix
# If per RG == 0 or undef, output is an actual file.
def basefreq(input, output, positions, per_rg):
    threads = config.param('bvatools_basefreq', 'threads', type='int')

    return Job(
        [input, positions],
        [output],
        [
            ['bvatools_basefreq', 'module_java'],
            ['bvatools_basefreq', 'module_bvatools']
        ],
        command="""\
java {java_other_options} -Xmx{ram} -jar $BVATOOLS_JAR \\
  basefreq \\
  --pos {positions} \\
  --bam {input}{per_rg}{threads} \\
  --out {output}""".format(
        java_other_options=config.param('bvatools_basefreq', 'java_other_options'),
        ram=config.param('bvatools_basefreq', 'ram'),
        positions=positions,
        per_rg=" \\\n  --per_rg " if per_rg else "",
        threads=" \\\n  --useIndex --threads " + str(threads) if threads > 1 else "",
        input=input,
        output=output
        ),
        removable_files=[output]
    )

def depth_of_coverage(input, output, coverage_bed, reference_genome="", other_options=""):
    return Job(
        [input, coverage_bed],
        [output],
        [
            ['bvatools_depth_of_coverage', 'module_java'],
            ['bvatools_depth_of_coverage', 'module_bvatools']
        ],
        command="""\
java {java_other_options} -Xmx{ram} -jar $BVATOOLS_JAR \\
  depthofcoverage {other_options} \\
  --threads {threads} \\
  --ref {reference_genome}{intervals} \\
  --bam {input} \\
  > {output}""".format(
        java_other_options=config.param('bvatools_depth_of_coverage', 'java_other_options'),
        ram=config.param('bvatools_depth_of_coverage', 'ram'),
        other_options=other_options,
        threads=config.param('bvatools_depth_of_coverage', 'threads', type='posint'),
        reference_genome=reference_genome if reference_genome else config.param('bvatools_depth_of_coverage', 'genome_fasta', type='filepath'),
        intervals=" \\\n  --intervals " + coverage_bed if coverage_bed else "",
        input=input,
        output=output
        )
    )

def extract_sclip(bamFile, output_prefix, flank="200"):
    return Job(
        [bamFile], 
        [
          output_prefix + ".sc.bam", 
          output_prefix + ".scOthers.bam", 
          output_prefix + ".scPositions.txt", 
          output_prefix + ".scSequences.txt"
        ], 
        [
          ['bvatools_ratiobaf', 'module_java'], 
          ['bvatools_ratiobaf', 'module_bvatools']
        ],
        command="""\
java {java_other_options} -Xmx{ram} -jar $BVATOOLS_JAR \\
  extractsclip {other_options} \\
  --bam {bamFile} \\
  --flank {flank} \\
  --minSCCount {minSCCount} \\
  --minSCLength {minSCLength} \\
  --minMappingQuality {minMappingQuality} \\
  --threads {threads} \\
  --prefix {output_prefix}""".format(
        java_other_options=config.param('extract_sclip', 'java_other_options'),
        ram=config.param('extract_sclip', 'ram'),
        other_options=config.param('extract_sclip', 'other_options', required=False),
        bamFile=bamFile,
        flank=flank,
        minSCCount=config.param('extract_sclip', 'min_sclip_count'),
        minSCLength=config.param('extract_sclip', 'kmer'),
        minMappingQuality=config.param('extract_sclip', 'min_mapping_quality'),
        threads=config.param('extract_sclip', 'threads'),
        output_prefix=output_prefix
        )
    )

def groupfixmate(input, output):
    return Job(
        [input],
        [output],
        [
            ['bvatools_groupfixmate', 'module_java'],
            ['bvatools_groupfixmate', 'module_bvatools']
        ],
        command="""\
java {java_other_options} -Xmx{ram} -jar $BVATOOLS_JAR \\
  groupfixmate \\
  --level 1 \\
  --bam {input} \\
  --out {output}""".format(
        java_other_options=config.param('bvatools_groupfixmate', 'java_other_options'),
        ram=config.param('bvatools_groupfixmate', 'ram'),
        input=input,
        output=output
        ),
        removable_files=[output]
    )

def ratiobaf(basefreq, output_prefix, positions):
    return Job(
        [basefreq, positions],
        [output_prefix + ".png"],
        [
            ['bvatools_ratiobaf', 'module_java'],
            ['bvatools_ratiobaf', 'module_bvatools']
        ],
        command="""\
java {java_other_options} -Xmx{ram} -jar $BVATOOLS_JAR \\
  ratiobaf {other_options} \\
  --refdict {reference_dictionary} \\
  --snppos {positions} \\
  --basefreq {basefreq} \\
  --prefix {output_prefix}""".format(
        java_other_options=config.param('bvatools_ratiobaf', 'java_other_options'),
        ram=config.param('bvatools_ratiobaf', 'ram'),
        other_options=config.param('bvatools_ratiobaf', 'other_options', required=False),
        reference_dictionary=config.param('bvatools_ratiobaf', 'genome_dictionary', type='filepath'),
        positions=positions,
        basefreq=basefreq,
        output_prefix=output_prefix
        )
    )

def readsqc(read1, read2, type, region_name, output_directory):
    threads = config.param('bvatools_readsqc', 'threads', type='int', required=False)

    return Job(
        [read1, read2],
        [output_directory + os.sep + "mpsQC_" + region_name + "_stats.xml"],
        [
            ['bvatools_readsqc', 'module_java'],
            ['bvatools_readsqc', 'module_bvatools']
        ],
        command="""\
java {java_other_options} -Xmx{ram} -jar $BVATOOLS_JAR \\
  readsqc {other_options} \\
  --regionName {region_name} \\
  --type {type} \\
  --output {output_directory} \\
  --read1 {read1}{read2}""".format(
        java_other_options=config.param('bvatools_readsqc', 'java_other_options'),
        ram=config.param('bvatools_readsqc', 'ram'),
        other_options=config.param('bvatools_readsqc', 'other_options', required=False),
        region_name=region_name,
        type=type,
        output_directory=output_directory,
        read1=read1,
        read2=" \\\n  --read2 " + read2 if read2 else "",
        threads=" \\\n  --threads " + str(threads) if threads > 1 else ""
        )
    )
