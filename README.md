# LoliumBlast
Methods and code for analyzing Lolium blast in Uruguay

## Ultra-accurate SNP calling
1. Generate SNP calls using Bowtie2/GATK
a. Align raw reads using bowtie2 (high sensitivity mode)
b. Call SNPs using GATK

2. Generate SNP calls using SNPcaller
a. Trim reads using Trimmomatic and assemble genomes using SPAdes
b. Standardize contig headers:
```
perl SimpleFastaHeaders.pl <GenomeID>.fasta
```
c. Trim short contigs and remove contaminating sequences using FCS:
```
genome=$1
prefix=${genome/_*/}
basename=${prefix/*\//}
run_fcsadaptor.sh --fasta-input $genome --output-dir $outdir \
 --prok --container-engine singularity --image $FCS_LOC/fcs-adaptor.sif
cat $genome | python3 /project/farman_s26abt480/public/fcs.py clean genome \
  --action-report ${basename}_FCS/fcs_adaptor_report.txt --output ${basename}_final.fasta
```
d. Mask repeat sequences using StrictRepeatMask:
```
perl StrictRepeatMask.pl $GenomeID_final.fasta
```
e. Align masked genome to masked reference:
```
blastn -query TF05-1_nh.fasta -subject <GenomeID>_final_masked.fasta -evalue 1e-20 -max_target_seqs 20000 -outft '6 qseqid sseqid qstart qend sstart send btop' > BLAST_DIR/TF05-1.<GenomeID>.BLAST
```
f. Call SNPs using SNPcaller:
```
perl SNPcaller.pl BLAST_DIR SNP_DIR
```
2. Validate SNP calls
a. Cross-reference GATK versus SNPcaller variants:
```
perl ValidateGatkSNPs.pl <GenomeID>_genotyped-snps.vcf SNP_DIR
```
b. Generate bed file listing validated SNP sitesa
```
perl SNPsummary.pl <GATK_DIR> | awk '{OFS="\t"; print $1, $2, $2, "SNP_" $2}' > All_SNP_sites.txt
```

