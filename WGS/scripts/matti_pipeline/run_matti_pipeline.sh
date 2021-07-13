# COPY FILES TO THE TAMBIO SERVER
rsync --progress -zrhv --bwlimit=5000000 --ignore-existing /groups/wyattgrp/data/bam/WGS_PCa/alignments/FILE_OF_INTEREST vpc@@tambio.uta.fi:DIR_OF_INTEREST/
#pass: ruby1867

# CREATE TUMOR-NORMAL PAIRS SHEET FOR THE COHORT
- The file must be tab-separated
- First line must contain column headers "TEST" and "REF"
- For tumor samples that have a matched control sample, add a line with the tumor sample in the first column, and the matched control in the second column
- For tumor samples without a matched control sample, add a line with the tumor sample in the first column, and an empty second column
- For control samples without any matched tumor samples, add a line with an empty first column, and the control sample in the second column

cd /home/vpc/datasets/mitochondrial_dna_project/alignments

# IDENTIFY SOMATIC MUTATIONS AND GERMLINE VARIANTS USING A TARGETED PANEL
echo X `seq 22` M | parallel -n8 'mutato call2 --region=chr${x} --alt-reads=4 --alt-frac=0.005 --max-frag-len=2000 ~/homo_sapiens/hg38.fa *.bam > ../mutations/chr${x}.vcf'
cat chr1.vcf <(cat chr{2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,M}.vcf | grep -v 'CHROM') > variants.vcf

# Somatic mutations
variant somatic --alt-reads=4 --alt-frac=0.005 --test-ref-ratio=3 --test-bg-ratio=20 --ref-reads=20 --min-sidedness=10 --min-mapq=10 variants.vcf ../tumor_normal_pairs.txt | variant predict effect - | variant protein altering - > somatic_protein_altering.tmp

variant nearby indels variants.vcf | variant somatic --alt-reads=10 --alt-frac=0.01 --test-ref-ratio=3 --test-bg-ratio=20 --ref-reads=20 --min-sidedness=30 --min-mapq=30 --mapq-filter-max-indel-len=100 - ../tumor_normal_pairs.txt | variant predict effect - | variant protein altering --invert - | variant discard sketchy silent - > somatic_silent.tmp

cat somatic_protein_altering.tmp <(tail -n +2 somatic_silent.tmp) | sort -k1,1V -k2,2n | variant annotate - ~/homo_sapiens/cosmic_77_hg38.jls | variant annotate - ~/homo_sapiens/exac_0.3_hg38.jls | variant discard if frequency above - ~/homo_sapiens/exac_0.3_hg38.jls 0.005 > somatic.vcf

# Germline variants
variant germline --alt-reads=8 --alt-frac=0.15 --bg-ratio=20 variants.vcf 'WBC|gDNA|enign' | variant predict effect - | variant protein altering - | variant annotate - ~/homo_sapiens/gnomad_v3.jls | variant discard if frequency above - ~/homo_sapiens/gnomad_v3.jls 0.005 | variant annotate - ~/homo_sapiens/cosmic_77_hg38.jls | variant annotate - ~/homo_sapiens/clinvar-2020-07-06_hg38.jls > germline.vcf

# Heterozygous germline SNPs
variant heterozygous snps --min-depth=30 variants.vcf 'WBC|gDNA|enign' | variant discard indels - | variant annotate - ~/homo_sapiens/gnomad_v3.jls | egrep 'CHROM|GNOMAD|ExAC' > hetz_snps.vcf
