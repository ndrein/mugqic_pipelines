#!/usr/bin/env bash

module load blat/35
module load mugqic-pipelines/2.2.0

blastdb=../reference-files
input_dir=bwa_search_with_contigs
output_dir=blat_search_with_contigs

mkdir $output_dir

blat -noHead -minIdentity=90 -minScore=50 $blastdb/microbial_all_cds.fasta $input_dir/cow_contigs_n_micro_cds.fasta -fine -q=rna -t=dna -out=blast8 $output_dir/cow_contigs.blatout

perl /hpf/largeprojects/ccmbio/nreinhardt/mugqic_pipelines/pipelines/metatranscriptomics/scripts/main_sort_blastout_fromfile.pl \
    $output_dir/cow_contigs_n_micro_cds.blatout \
    $output_dir/cow_conigs_n_micro_cds_sorted.blatout
#perl main_sort_blastout_fromfile.pl cow n_micro_cds blat contigs 10

perl /hpf/largeprojects/ccmbio/nreinhardt/mugqic_pipelines/pipelines/metatranscriptomics/scripts/main_get_blast_fromfile_1tophit.pl \
    $output_dir/cow_contigs_n_micro_cds_sorted.blatout \
    $output_dir/cow_contigs_n_micro_cds_blat_IDs.txt \
    $output_dir/cow_contigs_n_micro_cds_blat_pairs.txt \
    $output_dir/cow_contigs_n_micro_cds_blat_hitsID.txt \
    1 100 85 65 60
#perl main_get_blast_fromfile_1tophit.pl cow micro_cds blat contigs 1 100 85 65 60

python /hpf/largeprojects/ccmbio/nreinhardt/mugqic_pipelines/pipelines/metatranscriptomics/scripts/split_reads_by_id.py \
    --fastq $input_dir/cow_contigs_n_micro_cds.fasta \
    --id-file $output_dir/cow_contigs_n_micro_cds_blat_IDs.txt \
    --included $output_dir/cow_contigs_n_micro_cds_blat.fasta \
    --excluded $output_dir/cow_contigs_n_micro_cds_rest.fasta
#perl main_select_reads_fromfile.pl cow microgenes_blat blat contigs micro_cds

