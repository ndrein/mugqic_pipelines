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

        # output_prefix = 'format_reads'
        for readset in self.readsets:
            # output_dir = join(output_prefix, readset.name)
            #
            # output1 = join(output_dir, readset.name + '.1.formatted.fastq')
            # output2 = join(output_dir, readset.name + '.2.formatted.fastq')
            output1 = self.fm.declare_output('formatted-headers1', self.format_fastq_headers, {'readset_name': readset.name})
            output2 = self.fm.declare_output('formatted-headers2', self.format_fastq_headers, {'readset_name': readset.name})

            # sys.stderr.write(output1)
            # sys.exit()

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
                                    name='format_fastq_headers.' + readset.name))

        return jobs

    def trimmomatic(self):
        jobs = []

        # input_prefix = 'format_reads'
        # output_prefix = 'format_reads'
        #
        # def get_inputs(readset):
        #     """
        #     :return: 2 fastq filenames for paired-end reads
        #     """
        #     input_dir = join(input_prefix, readset.name)
        #     return input_dir, \
        #            join(input_dir, readset.name + '.1.formatted.fastq'), \
        #            join(input_dir, readset.name + '.2.formatted.fastq')
        #
        # def get_outputs(readset):
        #     """
        #     :return: output directory name,
        #              4 fastq filenames
        #     """
        #     output_dir = join(output_prefix, readset.name)
        #     return output_dir, \
        #            join(output_dir, readset.name + '.1.qual_paired.fastq'), \
        #            join(output_dir, readset.name + '.1.qual_unpaired.fastq'), \
        #            join(output_dir, readset.name + '.2.qual_paired.fastq'), \
        #            join(output_dir, readset.name + '.2.qual_unpaired.fastq')

        for readset in self.readsets:
            job = trimmomatic.trimmomatic(self.fm.find_output('formatted-headers1', self.format_fastq_headers),
                                          self.fm.find_output('formatted-headers2', self.format_fastq_headers),
                                          self.fm.declare_output('trim-paired1', self.trimmomatic, {'readset_name': readset.name}),
                                          self.fm.declare_output('trim-unpaired2', self.trimmomatic, {'readset_name': readset.name}),
                                          self.fm.declare_output('trim-paired2', self.trimmomatic, {'readset_name': readset.name}),
                                          self.fm.declare_output('trim-unpaired2', self.trimmomatic, {'readset_name': readset.name}),
                                          None,
                                          None,
                                          adapter_file=config.param('trimmomatic', 'adapter_fasta'),
                                          trim_log=join(self.trimmomatic.__name__, readset.name + '.trim.log'))
                                          # trim_log=join(output_prefix, readset.name + '.trim.log'))
            job.name = self.trimmomatic.__name__ + '.' + readset.name
            jobs.append(job)

        return jobs

    def merge_overlapping_reads(self):
        """
        Reads from the paired-end fastqs are merged together.

        The merged reads are added to the first paired-end fastq
        """
        jobs = []

        input_prefix = 'format_reads'
        output_prefix = 'format_reads'

        def get_inputs(readset):
            """
            :return: 2 fastq filenames for paired-end reads
            """
            input_dir = join(input_prefix, readset.name)
            return join(input_dir, readset.name + '.1.qual_paired.fastq'), \
                   join(input_dir, readset.name + '.2.qual_paired.fastq')

        def get_flash_params(readset):
            """
            :return: flash's output directory and output prefix
            """
            return readset.name, join(output_prefix, readset.name)

        def get_flash_outputs(output_dir, flash_output_prefix):
            """
            Get the filenames that flash will output

            Flash will output 3 files, one for both fastqs, and one for the merged reads
            :return: 2 uncombined fastq files + merged reads fastq file
            """
            return join(output_dir, flash_output_prefix + '.notCombined_1.fastq'), \
                   join(output_dir, flash_output_prefix + '.notCombined_2.fastq'), \
                   join(output_dir, flash_output_prefix + '.extendedFrags.fastq')

        def get_outputs(readset):
            """
            :return: 2 output filenames
            """
            output_dir = join(output_prefix, readset.name)
            return join(output_dir, readset.name + '.1.qual_all.fastq'), \
                   join(output_dir, readset.name + '.2.qual_all.fastq')

        for readset in self.readsets:
            # Get the filenames
            input1, input2 = get_inputs(readset)
            flash_output_dir, flash_output_prefix = get_flash_params(readset)
            flash1, flash2, flash_merged = get_flash_outputs(flash_output_dir, flash_output_prefix)
            output1, output2 = get_outputs(readset)

            # Create jobs
            flash_job = flash.merge_overlapping_reads(input1, input2, flash_output_dir, flash_output_prefix)
            # Put the merged reads into fastq1
            fastq1_job = Job(input_files=[flash1, flash_merged],
                             output_files=[output1],
                             command='cat {flash1} {merged} > {output1}'.format(flash1=flash1, merged=flash_merged,
                                                                                output1=output1))
            # Rename fastq2 to be consistent with fastq1
            fastq2_job = Job(input_files=[flash2],
                             output_files=[output2],
                             command='mv {flash2} {output2}'.format(flash2=flash2, output2=output2))
            jobs.append(concat_jobs([flash_job, fastq1_job, fastq2_job], name='flash.' + readset.name))

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

        input_prefix = 'format_reads'
        output_prefix = 'filter_reads'

        for readset in self.readsets:
            input_dir = join(input_prefix, readset.name)
            output_dir = join(output_prefix, readset.name)

            # For both fastq files (paired-end reads)
            for i in (1, 2):
                input_fastq = join(input_dir, readset.name + '.{i}.qual_all.fastq'.format(i=i))
                output_fasta = join(output_dir, readset.name + '.{i}.qual_all.fasta'.format(i=i))
                job_name = 'fastq_to_fasta.{name}.{i}'.format(name=readset.name, i=i)

                jobs.append(seqtk.fastq_to_fasta(input_fastq, output_fasta, name=job_name))

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
                input_unique_fasta = join(input_unique_dir, '{name}.{i}.usearch_out.fasta'.format(name=readset.name, i=i))
                input_uc = join(input_unique_dir, '{name}.{i}.usearch_out.uc'.format(name=readset.name, i=i))

                output_ids = join(output_dir, '{name}.{i}.cluster_sizes.json'.format(name=readset.name, i=i))
                output_fastq = join(output_dir, '{name}.{i}.unique.fastq'.format(name=readset.name, i=i))
                output_fasta = join(output_dir, '{name}.{i}.unique.fasta'.format(name=readset.name, i=i))

                job_name = 'remove_duplicates.{name}.{i}'.format(name=readset.name, i=i)

                jobs.append(Job(name=job_name,
                                input_files=[input_fastq, input_unique_fasta, input_uc],
                                output_files=[output_ids, output_fastq, output_fasta],
                                command='python {script_path}/remove_duplicates.py'.format(script_path=self.script_path)))

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
            self.fastq_to_fasta,
            self.cluster_duplicates,
            self.remove_duplicates, #6
            # self.remove_abundant_rrna,
        ]


if __name__ == '__main__':
    Metatranscriptomics()
