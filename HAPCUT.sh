spcodel="Rhisa1 Rhivi1 Suiame1 Suibr2 Suigr1"

for spcode in $spcodel
do
  ref=HAPCUT/${spcode}_MATA_denovo_concatenate_m.consensus.fasta
  python align_sam.py -i HD_mapping/${spcode}_MATA_denovo_concatenate.bam -r HDalignments/${spcode}_MATA_denovo_concatenate_ali.fasta -o HAPCUT/${spcode}_MATA_denovo_concatenate_m.sam
  samtools sort -n HAPCUT/${spcode}_MATA_denovo_concatenate_m.sam | samtools fixmate - - | samtools sort -o HAPCUT/${spcode}_MATA_denovo_concatenate_sm.bam
  bcftools mpileup -A -B -x -f $ref -d 1000 HAPCUT/${spcode}_MATA_denovo_concatenate_sm.bam | bcftools call -Ov -m > HAPCUT/${spcode}_MATA_HAPCUT.vcf
  extractHAIRS --bam HAPCUT/${spcode}_MATA_denovo_concatenate_sm.bam --VCF HAPCUT/${spcode}_MATA_HAPCUT.vcf --out HAPCUT/${spcode}_fragment_file
  HAPCUT2 --fragments HAPCUT/${spcode}_fragment_file --VCF HAPCUT/${spcode}_MATA_HAPCUT.vcf --output HAPCUT/${spcode}_HAPCUT
  bcftools view -Oz HAPCUT/${spcode}_HAPCUT.phased.VCF -o HAPCUT/${spcode}_HAPCUT.phased.vcf.gz
  bcftools index HAPCUT/${spcode}_HAPCUT.phased.vcf.gz
  rm HAPCUT/${spcode}_HapCut.fasta
  bcftools consensus -f $ref -H 1pIu HAPCUT/${spcode}_HAPCUT.phased.vcf.gz| \
  sed -E "s/^>.+$/>${spcode}_HapCut1/g" >> HAPCUT/${spcode}_HapCut.fasta
  bcftools consensus -f $ref -H 2pIu HAPCUT/${spcode}_HAPCUT.phased.vcf.gz| \
  sed -E "s/^>.+$/>${spcode}_HapCut2/g" >> HAPCUT/${spcode}_HapCut.fasta
  cat HDcontigs/${spcode}_MATA_denovo_concatenate.fasta >> HAPCUT/${spcode}_HapCut.fasta
  mafft HAPCUT/${spcode}_HapCut.fasta > HAPCUT/${spcode}_HapCut_compare.fasta
done

Rscript HAPCUT_conflict_test.R
