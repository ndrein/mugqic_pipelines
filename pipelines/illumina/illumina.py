#!/usr/bin/env python

# Python Standard Modules
import logging
import os
import re
import sys

# Append mugqic_pipeline directory to Python library path
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(sys.argv[0]))))

# MUGQIC Modules
from core.job import *
from core.pipeline import *
from bio.readset import *

from bio import picard
from bio import trimmomatic

log = logging.getLogger(__name__)

# Abstract pipeline gathering common features of all Illumina sequencing pipelines (readsets, trimming, etc.)
# Specific steps must be defined in Illumina children pipelines.
class Illumina(Pipeline):

    @property
    def readsets(self):
        return self._readsets

    @property
    def samples(self):
        return self._samples

    def sam_to_fastq(self):
        jobs = []
        for readset in self.readsets:
            if readset.bam and not readset.fastq1:
                if readset.run_type == "PAIRED_END":
                    readset.fastq1 = re.sub("\.bam$", ".pair1.fastq.gz", readset.bam)
                    readset.fastq2 = re.sub("\.bam$", ".pair2.fastq.gz", readset.bam)
                elif readset.run_type == "SINGLE_END":
                    readset.fastq1 = re.sub("\.bam$", ".single.fastq.gz", readset.bam)
                else:
                    raise Exception("Error: run type \"" + readset.run_type +
                    "\" is invalid for readset \"" + readset.name + "\" (should be PAIRED_END or SINGLE_END)!")

                job = picard.sam_to_fastq(readset.bam, readset.fastq1, readset.fastq2)
                job.name = "sam_to_fastq." + readset.name
                jobs.append(job)
        return jobs

    def trim(self):
        jobs = []
        for readset in self.readsets:
            trim_file_prefix = os.path.join("trim", readset.sample.name, readset.name + ".trim.")
            if readset.run_type == "PAIRED_END":
                if readset.bam and not readset.fastq1:
                    readset.fastq1 = re.sub("\.bam$", ".pair1.fastq.gz", readset.bam)
                    readset.fastq2 = re.sub("\.bam$", ".pair2.fastq.gz", readset.bam)
                job = trimmomatic.trimmomatic(
                    readset.fastq1,
                    readset.fastq2,
                    trim_file_prefix + "pair1.fastq.gz",
                    trim_file_prefix + "single1.fastq.gz",
                    trim_file_prefix + "pair2.fastq.gz",
                    trim_file_prefix + "single2.fastq.gz",
                    None,
                    readset.quality_offset,
                    trim_file_prefix + "out",
                    trim_file_prefix + "stats.csv"
                )
            elif readset.run_type == "SINGLE_END":
                if readset.bam and not readset.fastq1:
                    readset.fastq1 = re.sub("\.bam$", ".single.fastq.gz", readset.bam)
                job = trimmomatic.trimmomatic(
                    readset.fastq1,
                    None,
                    None,
                    None,
                    None,
                    None,
                    trim_file_prefix + "single.fastq.gz",
                    readset.quality_offset,
                    trim_file_prefix + "out",
                    trim_file_prefix + "stats.csv"
                )
            else:
                raise Exception("Error: run type \"" + readset.run_type +
                "\" is invalid for readset \"" + readset.name + "\" (should be PAIRED_END or SINGLE_END)!")
            job.name = "trim." + readset.name
            jobs.append(job)
        return jobs
