#!/bin/bash
# 12.05.16 Pieter Spealman
#!!!!!!!Define Set!!!!!!
# Set name will form the part of the namespace and directory structure.
	set_name=Carolina_03_02_23_SGD_CDS
	echo "Processing" ${set_name} "..."
#=============================
# 0. Load Modules
#=============================
	module load star/intel/2.7.6a
	module load samtools/intel/1.11
	module load cutadapt/3.1
	module load bedtools/intel/2.29.2
#=============================
# 1. Make directory names
#=============================
	echo "Building and setting directories..."
	work_dir=/scratch/ps163/STAR_${set_name}/
	data_dir=/scratch/cgsb/gresham/LABSHARE/Data/Ribosome_profiling/PS_Project_Carolino/fastq/
	fasta_dir=${work_dir}/fastq/
	star_idx_dir=${work_dir}/idx/
	rrna_idx=${star_idx_dir}rrna/
	genome_idx=${star_idx_dir}genome/
	tmp_dir=${work_dir}tmp/
	output_dir=${work_dir}output/
	QC_dir=${tmp_dir}QC/
	Processed_dir=${work_dir}/processed/
#=============================
# 2. Make directories
#=============================
	mkdir -p ${work_dir}
	mkdir -p ${data_dir}
	mkdir -p ${fasta_dir}
	mkdir -p ${star_idx_dir}
	mkdir -p ${rrna_idx}
	mkdir -p ${genome_idx}
	mkdir -p ${tmp_dir}
	mkdir -p ${output_dir}
	mkdir -p ${QC_dir}
	mkdir -p ${Processed_dir}
#=============================
# 3. Set annotations default parameters
#=============================
	echo "Setting Genome annotations ..."
	reference_dir=${genome_idx}
	cp /scratch/cgsb/gresham/pieter/genome/Saccharomyces_cerevisiae_SGD/S288C_reference_genome_R64-3-1_20210421/S288C_reference_sequence_R64-3-1_20210421.fsa ${reference_dir}
	genome_fa=${reference_dir}S288C_reference_sequence_R64-3-1_20210421.fsa
	#
	cp /scratch/cgsb/gresham/pieter/genome/Saccharomyces_cerevisiae_SGD/S288C_reference_genome_R64-3-1_20210421/saccharomyces_cerevisiae_R64-3-1_20210421.gff ${reference_dir}
	genome_gff=${reference_dir}saccharomyces_cerevisiae_R64-3-1_20210421.gff
	
	# cp /scratch/ps163/Carolina_03_18_2021/ensembl_50/Saccharomyces_cerevisiae.R64-1-1.ncrna.fa ${reference_dir}
	# rrna_fa=${reference_dir}Saccharomyces_cerevisiae.R64-1-1.ncrna.fa
	
	cp /scratch/cgsb/gresham/pieter/genome/Saccharomyces_cerevisiae_ensembl/Saccharomyces_cerevisiae.R64-1-1.ncrna_wo_ncrna_genes.fa ${reference_dir}
	rrna_fa=${reference_dir}/Saccharomyces_cerevisiae.R64-1-1.ncrna_wo_ncrna_genes.fa
	
#
		#==============================
		# 3a. RUN ONCE: Download Genome Annotations 
		#==============================
		# Run Once subsection. Downloads genome annotation files to the reference directory
		# echo "Downloading annotations for" ${set_name} "..."
		# mkdir -p ${reference_dir}
		# cd ${reference_dir}
		# rm -f *.fasta
		# rm -f *.gff
		# wget ###
		# mv *.fasta ${genome_fa}
		# wget ###
		# mv *.gff ${genome_gff}
		# cd ..
#
#==============================
# 4. Set Star Parameters
#==============================
# Sets parameters 
echo "Setting STAR parameters..."
	# set threads
	nproc=8
	# set number of mismatch
	nmismatch=10
	# set STAR alignments
	nc_align_params="--seedSearchLmax 10 --outFilterMultimapScoreRange 1 --outFilterMultimapNmax 100 --outFilterMismatchNmax ${nmismatch}  --alignIntronMax 100"
	align_params="--seedSearchLmax 10 --outFilterMultimapScoreRange 1 --outFilterScoreMin 20 --outFilterMultimapNmax 10 --outFilterMatchNmin 20 --outFilterMismatchNmax ${nmismatch} --outFilterIntronMotifs RemoveNoncanonical --scoreGap -8 --scoreGapNoncan -16 --alignIntronMax 100"
	# set output format
	SAM_params="--outSAMtype BAM Unsorted --outSAMmode NoQS --outSAMprimaryFlag AllBestScore"
#
		#==============================
		# 4a. OPTIONAL RUN: Test STAR parameters using genome annotations
		#==============================
		# Optional subsection. This tests STAR
		# ${STAR_bin}STAR --runThreadN 4 --runMode genomeGenerate --genomeDir ${genome_idx} --genomeFastaFiles ${genome_fa}
		# echo "If you got this far without an error, then it worked."
		# rm -Rf ${reference_dir}/data/
#
#==============================
# 5. Process fastq / SRA files
#==============================
# This section allows for the download of SRA files, or just direct access to local files.
# These files are then processed using cutadept to remove the adapter.
echo "Processing input sequence files ..."
#
# Adapter sequence: If this needs to be defined per sample, replace the ${adapter_seq} with sequence in each section
# AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
#adapter_seq_R1=CTGTAGGCACCATCAATATAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC
poly=AAAAAAAAAAAAA
adapter_seq_R1=${poly}AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA
	##_1
	input_file_a=RNA_DGY1657_A_R1
	output_file_1=RNA_DGY1657_R1
	temp_name=${output_file_1}
		cp ${data_dir}${input_file_a}.fastq.gz ${fasta_dir}${input_file_a}.fastq.gz
		gunzip -f ${fasta_dir}${input_file_a}.fastq.gz
		# #
		echo "Starting on " ${input_file_a} ${temp_name}
		cutadapt -a ${adapter_seq_R1} --cores=20 --revcomp -e 0.12 -m 12 -o ${fasta_dir}${temp_name}.fastq ${fasta_dir}${input_file_a}.fastq > ${fasta_dir}${temp_name}.fastq.log
		riboseq_fa_1=${fasta_dir}${temp_name}.fastq
		oriboprefix_1=${tmp_dir}${temp_name}_filter_
		griboprefix_1=${tmp_dir}${temp_name}_genome_
		transprefix_1=${tmp_dir}${temp_name}_transcript_
		riboseq_nrrna_fa_1=${oriboprefix_1}Unmapped.out.mate1
		#rm -f ${SRA_dir}${input_file_name}.fastq
# # #
	# ##_2
	input_file_a=RNA_DGY1657_B_R1
	output_file_2=RNA_DGY1657_R2
	#
	temp_name=${output_file_2}
		cp ${data_dir}${input_file_a}.fastq.gz ${fasta_dir}${input_file_a}.fastq.gz
		gunzip -f ${fasta_dir}${input_file_a}.fastq.gz
		#
		echo "Starting on " ${input_file_a} ${temp_name}
		cutadapt -a ${adapter_seq_R1} --cores=20 --revcomp -e 0.12 -m 12 -o ${fasta_dir}${temp_name}.fastq ${fasta_dir}${input_file_a}.fastq > ${fasta_dir}${temp_name}.fastq.log
		riboseq_fa_2=${fasta_dir}${temp_name}.fastq
		oriboprefix_2=${tmp_dir}${temp_name}_filter_
		griboprefix_2=${tmp_dir}${temp_name}_genome_
		transprefix_2=${tmp_dir}${temp_name}_transcript_
		riboseq_nrrna_fa_2=${oriboprefix_2}Unmapped.out.mate1
		# rm -f ${SRA_dir}${input_file_name}.fastq
		
# #SRR948563
	##_3
	input_file_a=RNA_DGY1728_C_R1
	output_file_3=RNA_DGY1728_R1
	#
	temp_name=${output_file_3}
		cp ${data_dir}${input_file_a}.fastq.gz ${fasta_dir}${input_file_a}.fastq.gz
		gunzip -f ${fasta_dir}${input_file_a}.fastq.gz
		
		echo "Starting on " ${input_file_a} ${temp_name}
		cutadapt -a ${adapter_seq_R1} --cores=20 --revcomp -e 0.12 -m 12 -o ${fasta_dir}${temp_name}.fastq ${fasta_dir}${input_file_a}.fastq > ${fasta_dir}${temp_name}.fastq.log
		riboseq_fa_3=${fasta_dir}${temp_name}.fastq
		oriboprefix_3=${tmp_dir}${temp_name}_filter_
		griboprefix_3=${tmp_dir}${temp_name}_genome_
		transprefix_3=${tmp_dir}${temp_name}_transcript_
		riboseq_nrrna_fa_3=${oriboprefix_3}Unmapped.out.mate1
		# rm -f ${SRA_dir}${input_file_name}.fastq 
		#
	# ##_4
	input_file_a=RNA_DGY1728_D_R1
	output_file_4=RNA_DGY1728_R2
	temp_name=${output_file_4}
		cp ${data_dir}${input_file_a}.fastq.gz ${fasta_dir}${input_file_a}.fastq.gz
		gunzip -f ${fasta_dir}${input_file_a}.fastq.gz
		
		echo "Starting on " ${input_file_a} ${temp_name}
		cutadapt -a ${adapter_seq_R1} --cores=20 --revcomp -e 0.12 -m 12 -o ${fasta_dir}${temp_name}.fastq ${fasta_dir}${input_file_a}.fastq > ${fasta_dir}${temp_name}.fastq.log
		riboseq_fa_4=${fasta_dir}${temp_name}.fastq
		oriboprefix_4=${tmp_dir}${temp_name}_filter_
		griboprefix_4=${tmp_dir}${temp_name}_genome_
		transprefix_4=${tmp_dir}${temp_name}_transcript_
		riboseq_nrrna_fa_4=${oriboprefix_4}Unmapped.out.mate1
		# rm -f ${SRA_dir}${input_file_name}.fastq 
		
	# # ##_5
	input_file_a=RNA_DGY1735_E_R1
	output_file_5=RNA_DGY1735_R1
	temp_name=${output_file_5}
		cp ${data_dir}${input_file_a}.fastq.gz ${fasta_dir}${input_file_a}.fastq.gz
		gunzip -f ${fasta_dir}${input_file_a}.fastq.gz
		
		echo "Starting on " ${input_file_a} ${temp_name}
		cutadapt -a ${adapter_seq_R1} --cores=20 --revcomp -e 0.12 -m 12 -o ${fasta_dir}${temp_name}.fastq ${fasta_dir}${input_file_a}.fastq > ${fasta_dir}${temp_name}.fastq.log
		riboseq_fa_5=${fasta_dir}${temp_name}.fastq
		oriboprefix_5=${tmp_dir}${temp_name}_filter_
		griboprefix_5=${tmp_dir}${temp_name}_genome_
		transprefix_5=${tmp_dir}${temp_name}_transcript_
		riboseq_nrrna_fa_5=${oriboprefix_5}Unmapped.out.mate1
		#rm -f ${SRA_dir}${input_file_name}.fastq
# #
	# ##_6
	input_file_a=RNA_DGY1735_H_R1
	output_file_6=RNA_DGY1735_R2
	temp_name=${output_file_6}
		cp ${data_dir}${input_file_a}.fastq.gz ${fasta_dir}${input_file_a}.fastq.gz
		gunzip -f ${fasta_dir}${input_file_a}.fastq.gz
		
		echo "Starting on " ${input_file_a} ${temp_name}
		cutadapt -a ${adapter_seq_R1} --cores=20 --revcomp -e 0.12 -m 12 -o ${fasta_dir}${temp_name}.fastq ${fasta_dir}${input_file_a}.fastq > ${fasta_dir}${temp_name}.fastq.log
		riboseq_fa_6=${fasta_dir}${temp_name}.fastq
		oriboprefix_6=${tmp_dir}${temp_name}_filter_
		griboprefix_6=${tmp_dir}${temp_name}_genome_
		transprefix_6=${tmp_dir}${temp_name}_transcript_
		riboseq_nrrna_fa_6=${oriboprefix_6}Unmapped.out.mate1
		# rm -f ${SRA_dir}${input_file_name}.fastq
		
#SRR948563
	##_7
	input_file_a=RNA_DGY_1726_J_R1
	output_file_7=RNA_DGY1726_R1
	temp_name=${output_file_7}
		cp ${data_dir}${input_file_a}.fastq.gz ${fasta_dir}${input_file_a}.fastq.gz
		gunzip -f ${fasta_dir}${input_file_a}.fastq.gz
		
		echo "Starting on " ${input_file_a} ${temp_name}
		cutadapt -a ${adapter_seq_R1} --cores=20 --revcomp -e 0.12 -m 12 -o ${fasta_dir}${temp_name}.fastq ${fasta_dir}${input_file_a}.fastq > ${fasta_dir}${temp_name}.fastq.log
		riboseq_fa_7=${fasta_dir}${temp_name}.fastq
		oriboprefix_7=${tmp_dir}${temp_name}_filter_
		griboprefix_7=${tmp_dir}${temp_name}_genome_
		transprefix_7=${tmp_dir}${temp_name}_transcript_
		riboseq_nrrna_fa_7=${oriboprefix_7}Unmapped.out.mate1
		# rm -f ${SRA_dir}${input_file_name}.fastq 
		#
	# ##_8
	input_file_a=RNA_DGY_1726_K_R1
	output_file_8=RNA_DGY1726_R2
	temp_name=${output_file_8}
		cp ${data_dir}${input_file_a}.fastq.gz ${fasta_dir}${input_file_a}.fastq.gz
		gunzip -f ${fasta_dir}${input_file_a}.fastq.gz
		
		echo "Starting on " ${input_file_a} ${temp_name}
		cutadapt -a ${adapter_seq_R1} --cores=20 --revcomp -e 0.12 -m 12 -o ${fasta_dir}${temp_name}.fastq ${fasta_dir}${input_file_a}.fastq > ${fasta_dir}${temp_name}.fastq.log
		riboseq_fa_8=${fasta_dir}${temp_name}.fastq
		oriboprefix_8=${tmp_dir}${temp_name}_filter_
		griboprefix_8=${tmp_dir}${temp_name}_genome_
		transprefix_8=${tmp_dir}${temp_name}_transcript_
		riboseq_nrrna_fa_8=${oriboprefix_8}Unmapped.out.mate1
		# rm -f ${SRA_dir}${input_file_name}.fastq
		
	# # ##_9
	input_file_a=RNA_DGY_1741_M_R1
	output_file_9=RNA_DGY1741_R1
	temp_name=${output_file_9}
		cp ${data_dir}${input_file_a}.fastq.gz ${fasta_dir}${input_file_a}.fastq.gz
		gunzip -f ${fasta_dir}${input_file_a}.fastq.gz
		
		echo "Starting on " ${input_file_a} ${temp_name}
		cutadapt -a ${adapter_seq_R1} --cores=20 --revcomp -e 0.12 -m 12 -o ${fasta_dir}${temp_name}.fastq ${fasta_dir}${input_file_a}.fastq > ${fasta_dir}${temp_name}.fastq.log
		riboseq_fa_9=${fasta_dir}${temp_name}.fastq
		oriboprefix_9=${tmp_dir}${temp_name}_filter_
		griboprefix_9=${tmp_dir}${temp_name}_genome_
		transprefix_9=${tmp_dir}${temp_name}_transcript_
		riboseq_nrrna_fa_9=${oriboprefix_9}Unmapped.out.mate1
		#rm -f ${SRA_dir}${input_file_name}.fastq
# #
	# ##_10
	input_file_a=RNA_DGY_1741_N_R1
	output_file_10=RNA_DGY1741_R2
	temp_name=${output_file_10}
		cp ${data_dir}${input_file_a}.fastq.gz ${fasta_dir}${input_file_a}.fastq.gz
		gunzip -f ${fasta_dir}${input_file_a}.fastq.gz
		
		echo "Starting on " ${input_file_a} ${temp_name}
		cutadapt -a ${adapter_seq_R1} --cores=20 --revcomp -e 0.12 -m 12 -o ${fasta_dir}${temp_name}.fastq ${fasta_dir}${input_file_a}.fastq > ${fasta_dir}${temp_name}.fastq.log
		riboseq_fa_10=${fasta_dir}${temp_name}.fastq
		oriboprefix_10=${tmp_dir}${temp_name}_filter_
		griboprefix_10=${tmp_dir}${temp_name}_genome_
		transprefix_10=${tmp_dir}${temp_name}_transcript_
		riboseq_nrrna_fa_10=${oriboprefix_10}Unmapped.out.mate1
		# rm -f ${SRA_dir}${input_file_name}.fastq
		
#SRR948563
	##_11
	input_file_a=RNA_DGY_1743_O_R1
	output_file_11=RNA_DGY1743_R1
	temp_name=${output_file_11}
		cp ${data_dir}${input_file_a}.fastq.gz ${fasta_dir}${input_file_a}.fastq.gz
		gunzip -f ${fasta_dir}${input_file_a}.fastq.gz
		
		echo "Starting on " ${input_file_a} ${temp_name}
		cutadapt -a ${adapter_seq_R1} --cores=20 --revcomp -e 0.12 -m 12 -o ${fasta_dir}${temp_name}.fastq ${fasta_dir}${input_file_a}.fastq > ${fasta_dir}${temp_name}.fastq.log
		riboseq_fa_11=${fasta_dir}${temp_name}.fastq
		oriboprefix_11=${tmp_dir}${temp_name}_filter_
		griboprefix_11=${tmp_dir}${temp_name}_genome_
		transprefix_11=${tmp_dir}${temp_name}_transcript_
		riboseq_nrrna_fa_11=${oriboprefix_11}Unmapped.out.mate1
		# rm -f ${SRA_dir}${input_file_name}.fastq 
		#
	# ##_12
	input_file_a=RNA_DGY_1743_P_R1
	output_file_12=RNA_DGY1743_R2
	temp_name=${output_file_12}
		cp ${data_dir}${input_file_a}.fastq.gz ${fasta_dir}${input_file_a}.fastq.gz
		gunzip -f ${fasta_dir}${input_file_a}.fastq.gz
		
		echo "Starting on " ${input_file_a} ${temp_name}
		cutadapt -a ${adapter_seq_R1} --cores=20 --revcomp -e 0.12 -m 12 -o ${fasta_dir}${temp_name}.fastq ${fasta_dir}${input_file_a}.fastq > ${fasta_dir}${temp_name}.fastq.log
		riboseq_fa_12=${fasta_dir}${temp_name}.fastq
		oriboprefix_12=${tmp_dir}${temp_name}_filter_
		griboprefix_12=${tmp_dir}${temp_name}_genome_
		transprefix_12=${tmp_dir}${temp_name}_transcript_
		riboseq_nrrna_fa_12=${oriboprefix_12}Unmapped.out.mate1
		# rm -f ${SRA_dir}${input_file_name}.fastq

#==============================
# 6. Begin STAR Run
#==============================
# This section loads the fastq files (Section 5) to STAR given the set parameters (Section 4)
# Load sequenced reads
	echo "Beginning" ${set_name} "STAR run ..."
		echo "Beginning" ${set_name} "STAR run ..."
	# #=============================
	# # step 6.b: align to transcriptome
	# #==============================
    echo "building rrna index..."
    STAR --runThreadN $nproc --runMode genomeGenerate --genomeDir ${rrna_idx} --genomeFastaFiles ${rrna_fa} --genomeSAindexNbases 5 --genomeChrBinNbits 11
	echo "filtering rRNA in RNA_seq..."
    STAR --runThreadN $nproc --genomeDir ${rrna_idx} --readFilesIn ${riboseq_fa_1} --outFileNamePrefix ${oriboprefix_1} --outStd SAM --outReadsUnmapped Fastx --outSAMmode NoQS ${nc_align_params} > /dev/null
    STAR --runThreadN $nproc --genomeDir ${rrna_idx} --readFilesIn ${riboseq_fa_2} --outFileNamePrefix ${oriboprefix_2} --outStd SAM --outReadsUnmapped Fastx --outSAMmode NoQS ${nc_align_params} > /dev/null
    STAR --runThreadN $nproc --genomeDir ${rrna_idx} --readFilesIn ${riboseq_fa_3} --outFileNamePrefix ${oriboprefix_3} --outStd SAM --outReadsUnmapped Fastx --outSAMmode NoQS ${nc_align_params} > /dev/null
    STAR --runThreadN $nproc --genomeDir ${rrna_idx} --readFilesIn ${riboseq_fa_4} --outFileNamePrefix ${oriboprefix_4} --outStd SAM --outReadsUnmapped Fastx --outSAMmode NoQS ${nc_align_params} > /dev/null
    STAR --runThreadN $nproc --genomeDir ${rrna_idx} --readFilesIn ${riboseq_fa_5} --outFileNamePrefix ${oriboprefix_5} --outStd SAM --outReadsUnmapped Fastx --outSAMmode NoQS ${nc_align_params} > /dev/null
    STAR --runThreadN $nproc --genomeDir ${rrna_idx} --readFilesIn ${riboseq_fa_6} --outFileNamePrefix ${oriboprefix_6} --outStd SAM --outReadsUnmapped Fastx --outSAMmode NoQS ${nc_align_params} > /dev/null
    STAR --runThreadN $nproc --genomeDir ${rrna_idx} --readFilesIn ${riboseq_fa_7} --outFileNamePrefix ${oriboprefix_7} --outStd SAM --outReadsUnmapped Fastx --outSAMmode NoQS ${nc_align_params} > /dev/null
    STAR --runThreadN $nproc --genomeDir ${rrna_idx} --readFilesIn ${riboseq_fa_8} --outFileNamePrefix ${oriboprefix_8} --outStd SAM --outReadsUnmapped Fastx --outSAMmode NoQS ${nc_align_params} > /dev/null
    STAR --runThreadN $nproc --genomeDir ${rrna_idx} --readFilesIn ${riboseq_fa_9} --outFileNamePrefix ${oriboprefix_9} --outStd SAM --outReadsUnmapped Fastx --outSAMmode NoQS ${nc_align_params} > /dev/null
    STAR --runThreadN $nproc --genomeDir ${rrna_idx} --readFilesIn ${riboseq_fa_10} --outFileNamePrefix ${oriboprefix_10} --outStd SAM --outReadsUnmapped Fastx --outSAMmode NoQS ${nc_align_params} > /dev/null
    STAR --runThreadN $nproc --genomeDir ${rrna_idx} --readFilesIn ${riboseq_fa_11} --outFileNamePrefix ${oriboprefix_11} --outStd SAM --outReadsUnmapped Fastx --outSAMmode NoQS ${nc_align_params} > /dev/null
    STAR --runThreadN $nproc --genomeDir ${rrna_idx} --readFilesIn ${riboseq_fa_12} --outFileNamePrefix ${oriboprefix_12} --outStd SAM --outReadsUnmapped Fastx --outSAMmode NoQS ${nc_align_params} > /dev/null
	#===
	# #=============================
	# # step 6.b: align to genome
	# #==============================
	echo "building genome index..."
	STAR --runThreadN $nproc --runMode genomeGenerate --genomeDir ${genome_idx} --genomeFastaFiles ${genome_fa} --genomeSAindexNbases 4
	echo "aligning ribo_seq to the genome..."
	STAR --runThreadN $nproc --genomeDir ${genome_idx} --readFilesIn ${riboseq_nrrna_fa_1} --outFileNamePrefix ${griboprefix_1} ${SAM_params} ${align_params}
	STAR --runThreadN $nproc --genomeDir ${genome_idx} --readFilesIn ${riboseq_nrrna_fa_2} --outFileNamePrefix ${griboprefix_2} ${SAM_params} ${align_params}
	STAR --runThreadN $nproc --genomeDir ${genome_idx} --readFilesIn ${riboseq_nrrna_fa_3} --outFileNamePrefix ${griboprefix_3} ${SAM_params} ${align_params}
	STAR --runThreadN $nproc --genomeDir ${genome_idx} --readFilesIn ${riboseq_nrrna_fa_4} --outFileNamePrefix ${griboprefix_4} ${SAM_params} ${align_params}
	STAR --runThreadN $nproc --genomeDir ${genome_idx} --readFilesIn ${riboseq_nrrna_fa_5} --outFileNamePrefix ${griboprefix_5} ${SAM_params} ${align_params}
	STAR --runThreadN $nproc --genomeDir ${genome_idx} --readFilesIn ${riboseq_nrrna_fa_6} --outFileNamePrefix ${griboprefix_6} ${SAM_params} ${align_params}
	STAR --runThreadN $nproc --genomeDir ${genome_idx} --readFilesIn ${riboseq_nrrna_fa_7} --outFileNamePrefix ${griboprefix_7} ${SAM_params} ${align_params}
	STAR --runThreadN $nproc --genomeDir ${genome_idx} --readFilesIn ${riboseq_nrrna_fa_8} --outFileNamePrefix ${griboprefix_8} ${SAM_params} ${align_params}
	STAR --runThreadN $nproc --genomeDir ${genome_idx} --readFilesIn ${riboseq_nrrna_fa_9} --outFileNamePrefix ${griboprefix_9} ${SAM_params} ${align_params}
	STAR --runThreadN $nproc --genomeDir ${genome_idx} --readFilesIn ${riboseq_nrrna_fa_10} --outFileNamePrefix ${griboprefix_10} ${SAM_params} ${align_params}
	STAR --runThreadN $nproc --genomeDir ${genome_idx} --readFilesIn ${riboseq_nrrna_fa_11} --outFileNamePrefix ${griboprefix_11} ${SAM_params} ${align_params}
	STAR --runThreadN $nproc --genomeDir ${genome_idx} --readFilesIn ${riboseq_nrrna_fa_12} --outFileNamePrefix ${griboprefix_12} ${SAM_params} ${align_params}

# ===
	echo "Completed" ${set_name} "STAR run ..."
#============
#	Clean up STAR
#=============
	# echo "Cleaning up" ${set_name} "STAR ..."
	# rm -f ${tmp_dir}*Unmapped*
	# rm -f ${tmp_dir}*_transcript*
	# rm -f ${tmp_dir}*_Aligned.toTranscriptome*
#=============
#	Quality Control
QC_suffix=_genome_Aligned.out
#=============
	echo "output located at ${output_dir}"
	echo "Beginning" ${set_name} "bam sorting and indexing ..."
	QC_name=${output_file_1}
	echo ${QC_name} "bam sorting and indexing ..."
	samtools sort -m 5000000000 ${tmp_dir}${QC_name}${QC_suffix}.bam -o ${output_dir}${QC_name}.sorted.bam
	samtools index ${output_dir}${QC_name}.sorted.bam
	samtools view -h ${output_dir}${QC_name}.sorted.bam -o ${output_dir}${QC_name}.sorted.sam
	rm -f ${tmp_dir}${QC_name}${QC_suffix}.bam
	bedtools coverage -counts -s -b ${output_dir}${QC_name}.sorted.bam -a ${genome_gff} > ${output_dir}${QC_name}.coverage.tab
	#===	_2
	QC_name=${output_file_2}
	echo ${QC_name} "bam sorting and indexing ..."
	samtools sort -m 5000000000 ${tmp_dir}${QC_name}${QC_suffix}.bam -o ${output_dir}${QC_name}.sorted.bam
	samtools index ${output_dir}${QC_name}.sorted.bam
	samtools view -h ${output_dir}${QC_name}.sorted.bam -o ${output_dir}${QC_name}.sorted.sam
	rm -f ${tmp_dir}${QC_name}${QC_suffix}.bam
	bedtools coverage -counts -s -b ${output_dir}${QC_name}.sorted.bam -a ${genome_gff} > ${output_dir}${QC_name}.coverage.tab
	#===	_3
	QC_name=${output_file_3}
	echo ${QC_name} "bam sorting and indexing ..."
	samtools sort -m 5000000000 ${tmp_dir}${QC_name}${QC_suffix}.bam -o ${output_dir}${QC_name}.sorted.bam
	samtools index ${output_dir}${QC_name}.sorted.bam
	samtools view -h ${output_dir}${QC_name}.sorted.bam -o ${output_dir}${QC_name}.sorted.sam
	rm -f ${tmp_dir}${QC_name}${QC_suffix}.bam
	bedtools coverage -counts -s -b ${output_dir}${QC_name}.sorted.bam -a ${genome_gff} > ${output_dir}${QC_name}.coverage.tab
	#===	_4
	QC_name=${output_file_4}
	echo ${QC_name} "bam sorting and indexing ..."
	samtools sort -m 5000000000 ${tmp_dir}${QC_name}${QC_suffix}.bam -o ${output_dir}${QC_name}.sorted.bam
	samtools index ${output_dir}${QC_name}.sorted.bam
	samtools view -h ${output_dir}${QC_name}.sorted.bam -o ${output_dir}${QC_name}.sorted.sam
	rm -f ${tmp_dir}${QC_name}${QC_suffix}.bam
	bedtools coverage -counts -s -b ${output_dir}${QC_name}.sorted.bam -a ${genome_gff} > ${output_dir}${QC_name}.coverage.tab
	# #
	echo "Beginning" ${set_name} "bam sorting and indexing ..."
	QC_name=${output_file_5}
	echo ${QC_name} "bam sorting and indexing ..."
	samtools sort -m 5000000000 ${tmp_dir}${QC_name}${QC_suffix}.bam -o ${output_dir}${QC_name}.sorted.bam
	samtools index ${output_dir}${QC_name}.sorted.bam
	samtools view -h ${output_dir}${QC_name}.sorted.bam -o ${output_dir}${QC_name}.sorted.sam
	rm -f ${tmp_dir}${QC_name}${QC_suffix}.bam
	bedtools coverage -counts -s -b ${output_dir}${QC_name}.sorted.bam -a ${genome_gff} > ${output_dir}${QC_name}.coverage.tab
	#===	_2
	QC_name=${output_file_6}
	echo ${QC_name} "bam sorting and indexing ..."
	samtools sort -m 5000000000 ${tmp_dir}${QC_name}${QC_suffix}.bam -o ${output_dir}${QC_name}.sorted.bam
	samtools index ${output_dir}${QC_name}.sorted.bam
	samtools view -h ${output_dir}${QC_name}.sorted.bam -o ${output_dir}${QC_name}.sorted.sam
	rm -f ${tmp_dir}${QC_name}${QC_suffix}.bam
	bedtools coverage -counts -s -b ${output_dir}${QC_name}.sorted.bam -a ${genome_gff} > ${output_dir}${QC_name}.coverage.tab
	# #===	_3
	QC_name=${output_file_7}
	echo ${QC_name} "bam sorting and indexing ..."
	samtools sort -m 5000000000 ${tmp_dir}${QC_name}${QC_suffix}.bam -o ${output_dir}${QC_name}.sorted.bam
	samtools index ${output_dir}${QC_name}.sorted.bam
	samtools view -h ${output_dir}${QC_name}.sorted.bam -o ${output_dir}${QC_name}.sorted.sam
	rm -f ${tmp_dir}${QC_name}${QC_suffix}.bam
	bedtools coverage -counts -s -b ${output_dir}${QC_name}.sorted.bam -a ${genome_gff} > ${output_dir}${QC_name}.coverage.tab
	#===	_4
	QC_name=${output_file_8}
	echo ${QC_name} "bam sorting and indexing ..."
	samtools sort -m 5000000000 ${tmp_dir}${QC_name}${QC_suffix}.bam -o ${output_dir}${QC_name}.sorted.bam
	samtools index ${output_dir}${QC_name}.sorted.bam
	samtools view -h ${output_dir}${QC_name}.sorted.bam -o ${output_dir}${QC_name}.sorted.sam
	rm -f ${tmp_dir}${QC_name}${QC_suffix}.bam
	bedtools coverage -counts -s -b ${output_dir}${QC_name}.sorted.bam -a ${genome_gff} > ${output_dir}${QC_name}.coverage.tab
	#===
	echo "Beginning" ${set_name} "bam sorting and indexing ..."
	QC_name=${output_file_9}
	echo ${QC_name} "bam sorting and indexing ..."
	samtools sort -m 5000000000 ${tmp_dir}${QC_name}${QC_suffix}.bam -o ${output_dir}${QC_name}.sorted.bam
	samtools index ${output_dir}${QC_name}.sorted.bam
	samtools view -h ${output_dir}${QC_name}.sorted.bam -o ${output_dir}${QC_name}.sorted.sam
	rm -f ${tmp_dir}${QC_name}${QC_suffix}.bam
	bedtools coverage -counts -s -b ${output_dir}${QC_name}.sorted.bam -a ${genome_gff} > ${output_dir}${QC_name}.coverage.tab
	#===	_2
	QC_name=${output_file_10}
	echo ${QC_name} "bam sorting and indexing ..."
	samtools sort -m 5000000000 ${tmp_dir}${QC_name}${QC_suffix}.bam -o ${output_dir}${QC_name}.sorted.bam
	samtools index ${output_dir}${QC_name}.sorted.bam
	samtools view -h ${output_dir}${QC_name}.sorted.bam -o ${output_dir}${QC_name}.sorted.sam
	rm -f ${tmp_dir}${QC_name}${QC_suffix}.bam
	bedtools coverage -counts -s -b ${output_dir}${QC_name}.sorted.bam -a ${genome_gff} > ${output_dir}${QC_name}.coverage.tab
	#===	_3
	QC_name=${output_file_11}
	echo ${QC_name} "bam sorting and indexing ..."
	samtools sort -m 5000000000 ${tmp_dir}${QC_name}${QC_suffix}.bam -o ${output_dir}${QC_name}.sorted.bam
	samtools index ${output_dir}${QC_name}.sorted.bam
	samtools view -h ${output_dir}${QC_name}.sorted.bam -o ${output_dir}${QC_name}.sorted.sam
	rm -f ${tmp_dir}${QC_name}${QC_suffix}.bam
	bedtools coverage -counts -s -b ${output_dir}${QC_name}.sorted.bam -a ${genome_gff} > ${output_dir}${QC_name}.coverage.tab
	#===	_4
	QC_name=${output_file_12}
	echo ${QC_name} "bam sorting and indexing ..."
	samtools sort -m 5000000000 ${tmp_dir}${QC_name}${QC_suffix}.bam -o ${output_dir}${QC_name}.sorted.bam
	samtools index ${output_dir}${QC_name}.sorted.bam
	samtools view -h ${output_dir}${QC_name}.sorted.bam -o ${output_dir}${QC_name}.sorted.sam
	rm -f ${tmp_dir}${QC_name}${QC_suffix}.bam
	bedtools coverage -counts -s -b ${output_dir}${QC_name}.sorted.bam -a ${genome_gff} > ${output_dir}${QC_name}.coverage.tab

#============
#	Clean up Set Run
#=============
echo "Cleaning up" ${set_name} "..."
	rm -rf ${tmp_dir}
	#mv ${tmp_dir}*.sorted.bam ${Processed_dir}
	#mv ${tmp_dir}*.bai ${Processed_dir}
	#rm -f ${fasta_dir}*.fastq
	#rm -f *.sam
	#rm -f *.bam
	#rm -f *.bedGraph
echo "Completed " ${set_name} "run."

