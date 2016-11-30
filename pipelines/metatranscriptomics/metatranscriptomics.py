#!/usr/bin/env python

# Python Standard Modules
import logging
import os.path
from os.path import join
import sys

# Append mugqic_pipelines directory to Python library path
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(sys.argv[0])))))

# MUGQIC Modules
from core.config import config
from core.job import Job, concat_jobs
from pipelines import common
from bfx import trimmomatic


log = logging.getLogger(__name__)


class Metatranscriptomics(common.Illumina):
    # Location for pipeline scripts, to be used by all steps
    # 'metatranscriptomics/scripts'
    script_path = os.path.join(os.path.dirname(__file__), 'scripts')

    def format_fastq_headers(self):
        input_dir = '../reference-files'
        input1 = input_dir + '/' + 'cow1.fastq'
        input2 = input_dir + '/' + 'cow2.fastq'

        output_dir = 'format_reads'
        output1 = output_dir + '/' + 'cow1_new.fastq'
        output2 = output_dir + '/' + 'cow2_new.fastq'

        return [concat_jobs([Job(command='mkdir {}'.format(output_dir)),
                             Job(name='format_fastq_headers',
                                 module_entries=[['DEFAULT', 'module_perl']],
                                 command='perl {script_path}/main_add_subID_reads_fastq.pl '
                                         '{input1} {output1} '
                                         '{input2} {output2}'.format(script_path=self.script_path,
                                                                     input1=input1,
                                                                     output1=output1,
                                                                     input2=input2,
                                                                     output2=output2))])]

    def trimmomatic(self):
        input_dir = 'format_reads'
        input1 = join(input_dir, 'cow1_new.fastq')
        input2 = join(input_dir, 'cow2_new.fastq')

        output_dir = 'format_reads'
        output_paired1 = join(output_dir, 'cow1_qual_paired.fastq')
        output_unpaired1 = join(output_dir, 'cow1_qual_unpaired.fastq')
        output_paired2 = join(output_dir, 'cow2_qual_paired.fastq')
        output_unpaired2 = join(output_dir, 'cow2_qual_unpaired.fastq')

        job = trimmomatic.trimmomatic(input1,
                                        input2,
                                        output_paired1,
                                        output_unpaired1,
                                        output_paired2,
                                        output_unpaired2,
                                        None,
                                        None,
                                        adapter_file=config.param('trimmomatic', 'adapter_fasta'),
                                        trim_log=join(output_dir, 'cow.trim.log'))
        job.name = 'trimmomatic.cow'
        return job

    def flash(self):
        return [concat_jobs([Job(command='mkdir -p ' + 'flash'),
                             Job(input_files=['trim/cow1_qual_paired.fastq',
                                              'trim/cow2_qual_paired.fastq'],
                                 output_files=['flash/cow_qual.extendedFrags.fastq',
                                               'flash/cow_qual.hist',
                                               'flash/cow_qual.histogram',
                                               'flash/cow_qual.notCombined_1.fastq',
                                               'flash/cow_qual.notCombined_2.fastq'],
                                 command='{flash} -M 75 -p 64 -t 2 -o cow_qual -d flash '
                                         'trim/cow1_qual_paired.fastq trim/cow2_qual_paired.fastq'.format(flash=config.param('flash', 'location'))),
                             Job(input_files=['flash/cow_qual.extendedFrags.fastq',
                                              'flash/cow_qual.notCombined_1.fastq'],
                                 output_files=['flash/cow1_qual_all.fastq'],
                                 command='cat flash/cow_qual.extendedFrags.fastq flash/cow_qual.notCombined_1.fastq > flash/cow1_qual_all.fastq'),
                             Job(input_files=['flash/cow_qual.notCombined_2.fastq'],
                                 output_files=['flash/cow2_qual_all.fastq'],
                                 command='cp flash/cow_qual.notCombined_2.fastq flash/cow2_qual_all.fastq')],
                            name='flash.cow')]

    def fastq_to_fasta(self):
        return [concat_jobs([Job(command='mkdir -p remove_duplicates'),
                             Job(input_files=['flash/cow1_qual_all.fastq',
                                              'flash/cow2_qual_all.fastq'],
                                 output_files=['remove_duplicates/cow1_qual_all.fasta',
                                               'remove_duplicates/cow2_qual_all.fasta'],
                                 module_entries=[
                                     ['fastq_to_fasta', 'module_seqtk'],
                                     ['fastq_to_fasta', 'module_perl']
                                 ],
                                 command='''\
seqtk seq -a flash/cow1_qual_all.fastq > remove_duplicates/cow1_qual_all.fasta
seqtk seq -a flash/cow2_qual_all.fastq > remove_duplicates/cow2_qual_all.fasta''')],
                            name='fastq_to_fasta.cow')]

    def remove_duplicates(self):
        return [
            Job(input_files=['remove_duplicates/cow1_qual_all.fasta',
                             'remove_duplicates/cow2_qual_all.fasta'],
                output_files=['remove_duplicates/cow1_qual_all_unique.fasta',
                              'remove_duplicates/cow1_qual_all_unique.uc',
                              'remove_duplicates/cow2_qual_all_unique.fasta',
                              'remove_duplicates/cow2_qual_all_unique.uc'],
                module_entries=[
                    ['remove_duplicates', 'module_usearch']
                ],
                command='''\
usearch --derep_fullseq --cluster remove_duplicates/cow1_qual_all.fasta --seedsout remove_duplicates/cow1_qual_all_unique.fasta --sizeout -uc remove_duplicates/cow1_qual_all_unique.uc
usearch --derep_fullseq --cluster remove_duplicates/cow2_qual_all.fasta --seedsout remove_duplicates/cow2_qual_all_unique.fasta --sizeout -uc remove_duplicates/cow2_qual_all_unique.uc
/hpf/largeprojects/ccmbio/nreinhardt/mugqic_pipelines/pipelines/metatranscriptomics/scripts/main_get_derepli_IDs.py''',
                name='remove_duplicates.cow')
        ]

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
                            'remove_duplicates/cow1_qual_all_unique.fasta'.format(rfam=config.param('remove_abundant_rrna', 'rfam_location'))),
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
                            'remove_duplicates/cow2_qual_all_unique.fasta'.format(rfam=config.param('remove_abundant_rrna', 'rfam_location')))
                ], name='remove_abundant_rrna.cow'
            )
        ]

    @property
    def steps(self):
        return [
            self.format_fastq_headers,
            self.trimmomatic,
            self.flash,
            self.fastq_to_fasta,
            self.remove_duplicates, #5
            self.remove_abundant_rrna,
        ]

if __name__ == '__main__':
    Metatranscriptomics()
