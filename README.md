# SUIMAT_pub

## align_sam.py
usage: python align_sam.py -i <input .sam file> -r <two aligned consensus sequences in .fasta> -o <output .sam file>
align and adjust position coordinates for two sam files
used in HAPCUT.sh
 
## HAPCUT_conflict_test.R
usage: Rscript HAPCUT_conflict_test.R
compare the results between HAPCUT2 and local de novo assemblies with the local de novo assembled sequences
used in HAPCUT.sh
 
## HAPCUT.sh
pipeline for phasing with HAPCUT2 and comparing the results with the local de novo assembled sequences
intermediate files in /HAPCUT
summarized to HAPCUT_summary.csv

## HAPCUT_summary.csv
summary of HAPCUT.sh

## /HAPCUT
intermediate files for HAPCUT.sh

## HDcontigs/
HD MAT alleles from different approaches
[species]_ref_concatenate.fasta: the sequences in reference genomes
[species]_denovo_concatenate.fasta: two alleles from locally de novo assembling
[species]_genome_concatenate.fasta: two alleles from mapping variant detection

## HDalignments/
[...]_ali.fasta: MAFFT alignment of the according sequences in HDcontigs/
HD1-HD2_denovo.fasta: collection of all locally de novo assembled sequences in species which are considered the haplotypes in the genomes 
HD1_aa.fasta: HD1 gene amino acid sequences from locally de novo assembled sequences in species
HD2_aa.fasta: HD2 gene amino acid sequences from locally de novo assembled sequences in species
HD1_CDS_aaali.fasta: HD1 CDS DNA alignment by condon (prank)
HD2_CDS_aaali.fasta: HD2 CDS DNA alignment by condon (prank)

## HD_mapping/
.bam file of reads mapped to according sequences in HDcontigs/


 
