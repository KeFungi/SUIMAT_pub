spcodel="Rhisa1 Rhivi1 Suiame1 Suibr2 Suigr1"

for spcode in $spcodel
do
  #aligned consensus sequences of two alleles
  ref=HAPCUT/${spcode}_MATA_denovo_concatenate_m.consensus.fasta
  
  #incorporate two .sam files by adjusting .sam coordinates
  python align_sam.py -i HD_mapping/${spcode}_MATA_denovo_concatenate.bam -r HDalignments/${spcode}_MATA_denovo_concatenate_ali.fasta -o HAPCUT/${spcode}_MATA_denovo_concatenate_m.sam
  
  #define variants between two alleles
  samtools sort -n HAPCUT/${spcode}_MATA_denovo_concatenate_m.sam | samtools fixmate - - | samtools sort -o HAPCUT/${spcode}_MATA_denovo_concatenate_sm.bam
  bcftools mpileup -A -B -x -f $ref -d 1000 HAPCUT/${spcode}_MATA_denovo_concatenate_sm.bam | bcftools call -Ov -m > HAPCUT/${spcode}_MATA_HAPCUT.vcf
  
  #phasing with HAPCUT2
  extractHAIRS --bam HAPCUT/${spcode}_MATA_denovo_concatenate_sm.bam --VCF HAPCUT/${spcode}_MATA_HAPCUT.vcf --out HAPCUT/${spcode}_fragment_file
  HAPCUT2 --fragments HAPCUT/${spcode}_fragment_file --VCF HAPCUT/${spcode}_MATA_HAPCUT.vcf --output HAPCUT/${spcode}_HAPCUT
  
  #extract phased sequences for comparison
  bcftools view -Oz HAPCUT/${spcode}_HAPCUT.phased.VCF -o HAPCUT/${spcode}_HAPCUT.phased.vcf.gz
  bcftools index HAPCUT/${spcode}_HAPCUT.phased.vcf.gz
  rm HAPCUT/${spcode}_HapCut.fasta
  bcftools consensus -f $ref -H 1pIu HAPCUT/${spcode}_HAPCUT.phased.vcf.gz| \
  sed -E "s/^>.+$/>${spcode}_HapCut1/g" >> HAPCUT/${spcode}_HapCut.fasta
  bcftools consensus -f $ref -H 2pIu HAPCUT/${spcode}_HAPCUT.phased.vcf.gz| \
  sed -E "s/^>.+$/>${spcode}_HapCut2/g" >> HAPCUT/${spcode}_HapCut.fasta
  cat HDcontigs/${spcode}_MATA_denovo_concatenate.fasta >> HAPCUT/${spcode}_HapCut.fasta
  
  #align two haplotypes from HAPCUT2 phasing and two haplotypes from local de novo assembling
  mafft HAPCUT/${spcode}_HapCut.fasta > HAPCUT/${spcode}_HapCut_compare.fasta
done

#count conflicts between HAPCUT2 phasing and local de novo assemblies
Rscript HAPCUT_conflict_test.R
