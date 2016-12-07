#!/usr/bin/env python
"""
Split fasta or fastq into 2 output files, based on a JSON file with IDs
"""
import argparse
import json
from Bio import SeqIO


# Parse args
def parse_args():
    arg_parser = argparse.ArgumentParser(description=__doc__)

    original_data = arg_parser.add_mutually_exclusive_group()
    original_data.add_argument('--fastq')
    original_data.add_argument('--fasta')

    arg_parser.add_argument('--id-file')

    arg_parser.add_argument('--included', help='Output filename for reads that are in the IDs')
    arg_parser.add_argument('--excluded', help='Output filename for reads that are not in the IDs')

    return arg_parser.parse_args()

def get_ids(id_file):
    """
    Get a set of ids from a JSON file

    Example JSON file:
    {
        "rows": [
                    {"id": "@SRR5"}
                    {"id": "@SRR7"}
                ]
    }

    Example output:
    {'@SRR5', '@SRR7'}

    :param id_file: JSON filename
    :return: set of str
    """
    with open(id_file) as f:
        return {row['id'] for row in json.load(f)['rows']}


def partition_reads(all_reads, ids):
    """
    Split all_reads into 2 sets based on where the read's id is in ids

    :param all_reads: iterable of Bio.SeqRecords
    :param ids: set of IDs
    :return: set, set
    """
    included, excluded = set(), set()

    for read in all_reads:
        included.add(read) if read.id in ids else excluded.add(read)

    return included, excluded


def get_all_reads(file, format):
    return set(SeqIO.parse(file, format))


args = parse_args()
input, format = (args.fastq, 'fastq') if args.fastq else (args.fasta, 'fasta')

included, excluded = partition_reads(get_all_reads(input, format), get_ids(args.id_file))

SeqIO.write(included, args.included, format)
SeqIO.write(excluded, args.excluded, format)