#!/usr/bin/env python

################################################################################
# Copyright (C) 2014, 2015 GenAP, McGill University and Genome Quebec Innovation Centre
#
# This file is part of MUGQIC Pipelines.
#
# MUGQIC Pipelines is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# MUGQIC Pipelines is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with MUGQIC Pipelines.  If not, see <http://www.gnu.org/licenses/>.
################################################################################

# Python Standard Modules
import os

# MUGQIC Modules
from core.config import *
from core.job import *

def deseq(
    design_file,
    count_matrix,
    output_dir
    ):

    return  Job(
        [count_matrix],
        [os.path.join(output_dir, "deseq_results.csv"), os.path.join(output_dir, "dge_results.csv")],
        [
            ['differential_expression_deseq', 'module_mugqic_tools'],
            ['differential_expression_deseq', 'module_R']
        ],
        command="""\
Rscript $R_TOOLS/deseq.R \\
  -d {design_file} \\
  -c {count_matrix} \\
  -o {output_dir}""".format(
        design_file=design_file,
        count_matrix=count_matrix,
        output_dir=output_dir
    ))

def edger(
    design_file,
    count_matrix,
    output_dir
    ):

    return  Job(
        [count_matrix],
        [os.path.join(output_dir, "edger_results.csv")],
        [
            ['differential_expression_edger', 'module_mugqic_tools'],
            ['differential_expression_edger', 'module_R']
        ],
        command="""\
Rscript $R_TOOLS/edger.R \\
  -d {design_file} \\
  -c {count_matrix} \\
  -o {output_dir}""".format(
        design_file=design_file,
        count_matrix=count_matrix,
        output_dir=output_dir
    ))

def goseq(
    input_file,
    input_columns,
    output_file,
    ):

    return  Job(
        [input_file],
        [output_file],
        [
            ['differential_expression_goseq', 'module_mugqic_tools'],
            ['differential_expression_goseq', 'module_R']
        ],
        command="""\
Rscript $R_TOOLS/goseq.R {other_options} \\
  -a {gene_size_file} \\
  -G {gene_ontology_file} \\
  -d {input_file} \\
  -c {input_columns} \\
  -o {output_file}""".format(
        other_options=config.param('differential_expression_goseq','other_options'),
        gene_size_file=config.param('differential_expression_goseq', 'gene_size', type='filepath'),
        gene_ontology_file=config.param('differential_expression_goseq', 'gene_ontology', type='filepath'),
        input_file=input_file,
        input_columns=input_columns,
        output_file=output_file
    ))
