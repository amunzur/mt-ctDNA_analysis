#!/bin/bash

#SBATCH --job-name=Alignment-mito
#SBATCH -p debug,express,normal,big-mem,long
#SBATCH --cpus-per-task=5
#SBATCH --mem 200000 # memory pool for all cores
#SBATCH -t 4:00:00 # time (D-HH:MM or HH:MM:SS)
#SBATCH --export=all
#SBATCH --output=/groups/wyattgrp/log/%j.log
#SBATCH --error=/groups/wyattgrp/log/%j.log

bam_file=$1
main_mito_path="/groups/wyattgrp/users/amunzur/mt-ctDNA/WGS"
R1_fasta="${main_mito_path}/processed_fasta/${bam_file//.bam/.R1.fasta}" # string replacement
R2_fasta=${R1_fasta/R1/R2}

# NUCLEAR GENOME PARAMS
HG_sam="${main_mito_path}/sams/realigned_HG/${bam_file//.bam/.sam}" # sam file generated after aligning to human genome 
HG_bam="${main_mito_path}/bams/realigned_HG/${bam_file}" 
HG_bam_sorted="${HG_bam//.bam/_sorted.bam}" 
path_to_ref_HG="/groups/wyattgrp/users/amunzur/chip_project/references/hg38_no_mito.fa" # chrM removed from hg38
path_to_flagstat_HG="${main_mito_path}/metrics/flagstat/nuclear_genome_alignment/${bam_file//.bam/.txt}"

# MITO GENOME PARAMS
mito_sam="${main_mito_path}/sams/realigned_mito/${bam_file//.bam/.sam}" # sam file generated after aligning to human genome 
mito_bam="${main_mito_path}/bams/realigned_mito/${bam_file}"  
mito_bam_sorted="${mito_bam//.bam/_sorted.bam}" 
path_to_ref_mito="/groups/wyattgrp/users/amunzur/mt-ctDNA/references/hg38_mito.fa"
path_to_flagstat_mito="${main_mito_path}/metrics/flagstat/mito_genome_alignment/${bam_file//.bam/.txt}"

# bowtie stuff
alignment_tool="bowtie2"
bowtie_idx_path_nuclear="${path_to_ref_HG//hg38_no_mito.fa/hg38}"
bowtie_idx_path_mito="/groups/wyattgrp/users/amunzur/mt-ctDNA/bowtie_idx/hg38_mito_bowtie.fa"
bowtie_log_main="${main_mito_path}/log/bowtie2/" # warnings etc will go here as a text file

# align and index
align_and_samtools () {

	# both alignment tools output a sam file
	if [ $8 = "bwa" ]
	then 
		echo "Started aligning with BWA."

		bwa mem -t 6 -R "@RG\tID:sample" $1 $2 $3 > $4

		samtools view -bS $4 > $5 # sam to bam
		samtools sort $5 -o $6 # sort the bam
		samtools index $6 # index the sorted bam
		samtools flagstat $6 > $7 # run flagstat metrics \

	else
		echo "Started aligning with Bowtie2."

		# doesnt matter if this is nuclear or mito
		sam="${4//sams/sams/bowtie_sams}" # sam file generated after aligning to human genome 
		bam="${5//bams/bams/bowtie_bams}" # bam file converted from sam
		sorted_bam="${bam//.bam/_sorted.bam}"
		path_to_flagstat="${7//metrics/bowtie_metrics}"

		bowtie2 -f -X 1000 -x $9 -1 $2 -2 $3 1> $sam 2> $bowtie_log

		samtools view -bS $sam > $bam # sam to bam
		samtools sort $bam -o $sorted_bam # sort the bam
		samtools index $6 # index the sorted bam
		samtools flagstat $6 > $path_to_flagstat # run flagstat metrics 

	fi
	echo "Finished aligning."

}

# echo "ALIGNING TO THE NUCLEAR GENOME."
# echo  $bowtie_idx_path_nuclear
# align_and_samtools $path_to_ref_HG $R1_fasta $R2_fasta $HG_sam $HG_bam $HG_bam_sorted $path_to_flagstat_HG $alignment_tool $bowtie_idx_path_nuclear $bowtie_log_main# nuclear

echo "ALIGNING TO THE MITO GENOME."
bowtie_log="${bowtie_log_main}${bam_file//.bam/.txt}"
align_and_samtools $path_to_ref_mito $R1_fasta $R2_fasta $mito_sam $mito_bam $mito_bam_sorted $path_to_flagstat_mito $alignment_tool $bowtie_idx_path_mito $bowtie_log # mito

# echo "PROCESSING THE NUCLEAR ALIGNMENT."
# filter_with_GATK $HG_bam filtered_${HG_bam}

# echo "PROCESSING THE MITO ALIGNMENT."





