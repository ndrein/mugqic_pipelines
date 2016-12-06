#!/usr/bin/env python

# Python Standard Modules
import logging
import os.path
import sys
from os.path import join

# Append mugqic_pipelines directory to Python library path
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(sys.argv[0])))))

# MUGQIC Modules
from core.config import config
from core.job import Job, concat_jobs, mkdir
from pipelines import common
from bfx import trimmomatic
from bfx import flash
from bfx import seqtk
from bfx import usearch
from core.filename_manager import FilenameManager

log = logging.getLogger(__name__)


class Metatranscriptomics(common.Illumina):
    """
    Based off of the pipeline located here:
    http://bioinformatics-ca.github.io/analysis_of_metagenomic_data_module6_lab_2016/

    There are 4 phases to this pipeline so far:
    * format_reads
        Format read headers,
        trim reads,
        merge reads together
    * filter_reads
        Remove duplicates (temporarily),
        remove rRNA,
        remove reads from host
    * contig_assembly
        Assemble reads into larger contigs,
        map the reads to these contigs
    * search
        Search known microbial sequences using both contigs and singleton reads,
        try bwa, blat, and diamond
    * Further
    Each phase is associated with its own directory
    Within each phase directory, each sample/readset will have its own directory
    """

    # Location for pipeline scripts, to be used by all steps
    # 'metatranscriptomics/scripts'
    script_path = os.path.join(os.path.dirname(__file__), 'scripts')

    def __init__(self):
        self.fm = FilenameManager(self.steps)
        super(Metatranscriptomics, self).__init__()

    def format_fastq_headers(self):
        """
        Mark the headers of the fastq files with /1 or /2 to differentiate the paired-end reads

        Call 'main_add_subID_reads_fastq.pl'
        """
        jobs = []

        for readset in self.readsets:
            output1 = self.fm.declare_output('formatted-headers1.fastq', self.format_fastq_headers, readset.name)
            output2 = self.fm.declare_output('formatted-headers2.fastq', self.format_fastq_headers, readset.name)

            jobs.append(concat_jobs([mkdir(output1),
                                     mkdir(output2),
                                     Job(input_files=[readset.fastq1, readset.fastq2],
                                         output_files=[output1, output2],
                                         module_entries=[['DEFAULT', 'module_perl']],
                                         command='perl {script_path}/main_add_subID_reads_fastq.pl '
                                                 '{input1} {output1} '
                                                 '{input2} {output2}'.format(script_path=self.script_path,
                                                                             input1=readset.fastq1,
                                                                             output1=output1,
                                                                             input2=readset.fastq2,
                                                                             output2=output2))],
                                    name='{step}.{readset}'.format(step=self.format_fastq_headers.__name__,
                                                                   readset=readset.name)))

        return jobs

    def trimmomatic(self):
        jobs = []

        for readset in self.readsets:
            jobs.append(concat_jobs([trimmomatic.trimmomatic(
                self.fm.find_output('formatted-headers1.fastq', self.format_fastq_headers, readset.name),
                self.fm.find_output('formatted-headers2.fastq', self.format_fastq_headers, readset.name),
                self.fm.declare_output('trim-paired1.fastq', self.trimmomatic, readset.name),
                self.fm.declare_output('trim-unpaired2.fastq', self.trimmomatic, readset.name),
                self.fm.declare_output('trim-paired2.fastq', self.trimmomatic, readset.name),
                self.fm.declare_output('trim-unpaired2.fastq', self.trimmomatic, readset.name),
                None,
                None,
                adapter_file=config.param('trimmomatic', 'adapter_fasta'),
                trim_log=join(self.trimmomatic.__name__, readset.name + '.trim.log')),
                Job(command='mkdir -p {}'.format(self.trimmomatic.__name__))],
                name='{step}.{readset}'.format(step=self.trimmomatic.__name__, readset=readset.name)
            ))

        return jobs

    def merge_overlapping_reads(self):
        """
        Reads from the paired-end fastqs are merged together.

        """
        jobs = []

        flash_output_prefix = 'out'

        for readset in self.readsets:
            output_dir = join(self.merge_overlapping_reads.__name__, readset.name)

            jobs.append(flash.merge_overlapping_reads(
                fastq1=self.fm.find_output('trim-paired1.fastq', self.trimmomatic, readset.name),
                fastq2=self.fm.find_output('trim-paired2.fastq', self.trimmomatic, readset.name),
                output_dir=output_dir, output_prefix=flash_output_prefix))
            self.fm.declare_output(flash_output_prefix + '.notCombined_1.fastq',
                                   self.merge_overlapping_reads,
                                   readset.name)
            self.fm.declare_output(flash_output_prefix + '.notCombined_2.fastq',
                                   self.merge_overlapping_reads,
                                   readset.name)
            self.fm.declare_output(flash_output_prefix + '.extendedFrags.fastq',
                                   self.merge_overlapping_reads,
                                   readset.name)

        return jobs

    def merge_flash_files(self):
        """
        The merged reads are added to the first paired-end fastq
        """
        jobs = []

        for readset in self.readsets:
            flash_output_prefix = readset.name

            # Put the merged reads into fastq1
            flash1 = self.fm.find_output(flash_output_prefix + '.notCombined_1.fastq',
                                         self.merge_overlapping_reads,
                                         readset.name)
            flash_merged = self.fm.find_output(flash_output_prefix + '.extendedFrags.fastq',
                                               self.merge_overlapping_reads,
                                               readset.name)
            output1 = self.fm.declare_output('flash-1.fastq', self.merge_overlapping_reads, readset.name)
            fastq1_job = Job(input_files=[flash1,
                                          flash_merged],
                             output_files=[output1],
                             command='cat {flash1} {merged} > {output1}'.format(flash1=flash1,
                                                                                merged=flash_merged,
                                                                                output1=output1))

            # Rename fastq2 to be consistent with fastq1
            flash2 = self.fm.find_output(flash_output_prefix + '.notCombined_2.fastq',
                                         self.merge_overlapping_reads,
                                         readset.name)
            output2 = self.fm.declare_output('flash-2.fastq', self.merge_overlapping_reads, readset.name)
            fastq2_job = Job(input_files=[flash2],
                             output_files=[output2],
                             command='mv {flash2} {output2}'.format(flash2=flash2, output2=output2))
            jobs.append(concat_jobs([fastq1_job, fastq2_job], name='{step}.{name}'.format(step=self.merge_flash_files.__name__, name=readset.name)))

        return jobs

    def fastq_to_fasta(self):
        """
        Convert both fastq files to fastas

        Input:
        format_reads/*{1,2}.qual_all.fastq
        Output:
        filter_reads/*{1,2}.qual_all.fasta

        Required since our usearch/5.2.236 only takes fastas as input
        """
        jobs = []

        # input_prefix = 'format_reads'
        # output_prefix = 'filter_reads'

        for readset in self.readsets:
            # input_dir = join(input_prefix, readset.name)
            # output_dir = join(output_prefix, readset.name)

            # For both fastq files (paired-end reads)
            for i in (1, 2):
                # input_fastq = join(input_dir, readset.name + '.{i}.qual_all.fastq'.format(i=i))
                # output_fasta = join(output_dir, readset.name + '.{i}.qual_all.fasta'.format(i=i))
                # job_name = 'fastq_to_fasta.{name}.{i}'.format(name=readset.name, i=i)

                jobs.append(seqtk.fastq_to_fasta(self.fm.find_output('flash-{}.fastq'.format(i),
                                                                     self.merge_overlapping_reads,
                                                                     readset.name),
                                                 self.fm.declare_output('flash-{}.fasta'.format(i),
                                                                        self.fastq_to_fasta,
                                                                        readset.name),
                                                 name='{step}.{name}'.format(step=self.fastq_to_fasta.__name__,
                                                                             name=readset.name)))
                # jobs.append(seqtk.fastq_to_fasta(input_fastq, output_fasta, name=job_name))

        return jobs

    def cluster_duplicates(self):
        """
        Cluster duplicate reads together

        Input:
        filter_reads/*.{1,2}.qual_all.fasta

        Output:
        filter_reads/*.{1,2}.qual_all_unique.fasta
                     *.{1,2}.qual_all_unique.uc         - contains the fasta IDs and cluster IDs
        """
        jobs = []

        input_prefix = 'filter_reads'
        output_prefix = 'filter_reads'

        for readset in self.readsets:
            input_dir = join(input_prefix, readset.name)
            output_dir = join(output_prefix, readset.name)

            # For both paired-end reads
            for i in (1, 2):
                input_fasta = join(input_dir, '{name}.{i}.qual_all.fasta'.format(name=readset.name, i=i))

                output_fasta = join(output_dir, '{name}.{i}.usearch_out.fasta'.format(name=readset.name, i=i))
                output_uc = join(output_dir, '{name}.{i}.usearch_out.uc'.format(name=readset.name, i=i))

                job_name = 'cluster_duplicates.{name}.{i}'.format(name=readset.name, i=i)

                jobs.append(usearch.cluster_duplicates(input_fasta, output_fasta, output_uc, name=job_name))

        return jobs

    def remove_duplicates(self):
        jobs = []

        input_original_prefix = 'format_reads'
        input_unique_prefix = 'filter_reads'
        output_prefix = 'filter_reads'

        for readset in self.readsets:
            input_original_dir = join(input_original_prefix, readset.name)
            input_unique_dir = join(input_unique_prefix, readset.name)
            output_dir = join(output_prefix, readset.name)

            for i in (1, 2):
                input_fastq = join(input_original_dir, '{name}.{i}.qual_all.fastq'.format(name=readset.name, i=i))
                input_unique_fasta = join(input_unique_dir,
                                          '{name}.{i}.usearch_out.fasta'.format(name=readset.name, i=i))
                input_uc = join(input_unique_dir, '{name}.{i}.usearch_out.uc'.format(name=readset.name, i=i))

                output_ids = join(output_dir, '{name}.{i}.cluster_sizes.json'.format(name=readset.name, i=i))
                output_fastq = join(output_dir, '{name}.{i}.unique.fastq'.format(name=readset.name, i=i))
                output_fasta = join(output_dir, '{name}.{i}.unique.fasta'.format(name=readset.name, i=i))

                job_name = 'remove_duplicates.{name}.{i}'.format(name=readset.name, i=i)

                jobs.append(Job(name=job_name,
                                input_files=[input_fastq, input_unique_fasta, input_uc],
                                output_files=[output_ids, output_fastq, output_fasta],
                                command='python {script_path}/remove_duplicates.py'.format(
                                    script_path=self.script_path)))

        return jobs

    def remove_abundant_rrna(self):
        return [
            concat_jobs([
                Job(command='mkdir -p remove_abundant_rrna'),
                Job(input_files=['remove_duplicates/cow1_qual_all_unique.fasta'],
                    output_files=['remove_abundant_rrna/cow1_rRNA.log',
                                  'remove_abundant_rrna/cow1_rRNa.infernalout'],
                    module_entries=[
                        ['remove_abundant_rrna', 'module_infernal'],
                        ['remove_abundant_rrna', 'module_perl']
                    ],
                    command='cmscan -o remove_abundant_rrna/cow1_rRNA.log '
                            '--tblout remove_abundant_rrna/cow1_rRNA.infernalout '
                            '--noali --notextw --rfam -E 0.001 '
                            '{rfam} '
                            'remove_duplicates/cow1_qual_all_unique.fasta'.format(
                        rfam=config.param('remove_abundant_rrna', 'rfam_location'))),
                Job(input_files=['remove_duplicates/cow2_qual_all_unique.fasta'],
                    output_files=['remove_abundant_rrna/cow2_rRNA.log',
                                  'remove_abundant_rrna/cow2_rRNa.infernalout'],
                    module_entries=[
                        ['remove_abundant_rrna', 'module_infernal'],
                        ['remove_abundant_rrna', 'module_perl']
                    ],
                    command='cmscan -o remove_abundant_rrna/cow2_rRNA.log '
                            '--tblout remove_abundant_rrna/cow2_rRNA.infernalout '
                            '--noali --notextw --rfam -E 0.001 '
                            '{rfam} '
                            'remove_duplicates/cow2_qual_all_unique.fasta'.format(
                        rfam=config.param('remove_abundant_rrna', 'rfam_location')))
            ], name='remove_abundant_rrna.cow'
            )
        ]

    @property
    def steps(self):
        return [
            self.format_fastq_headers,
            self.trimmomatic,
            self.merge_overlapping_reads,  # 3
            self.merge_flash_files,
            self.fastq_to_fasta,
            self.cluster_duplicates,
            self.remove_duplicates,  # 7
            # self.remove_abundant_rrna,
        ]


if __name__ == '__main__':
    Metatranscriptomics()
