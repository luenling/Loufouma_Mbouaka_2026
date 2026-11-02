#!/bin/bash
#----------
# author: Lukas Endler
# authored: Friday, May 24, 2024 @ 12:27:00 PM
# Time-stamp: <Friday, May 24, 2024 @ 12:33:06 PM lukasendler>
# command lines for cufflinks
#--------------

WORKDIR=$(pwd)
echo WORKDIR=$WORKDIR
SCRIPTDIR=$(dirname -- "$0")
echo SCRIPTDIR=$SCRIPTDIR

BASEDIR=/shared_data/shared_data_BACKUP/experimental_data/P_2024_ACRS_AcanthamoebaRNASeq_LE
RAWDATA=${BASEDIR}/Reads/
REFBASE=${BASEDIR}/References/
BT2IDX="/shared_data/shared_data_BACKUP/auxiliary_data/indices/Human/GRCh38.p14/GRCh38_noalt_as_bowtie2/GRCh38_noalt_as"
SCRIPTS=/${BASEDIR}/Scripts/
LOGFILE=${WORKDIR}/"logfile.log"
ERRLOG=${WORKDIR}/"logfile.err.log"

source "${SCRIPTDIR}/bash_functions.sh"

mkdir -p ${WORKDIR}/LOGS
mkdir -p ${WORKDIR}/Reports
BASENAMES=$@

declare -A GENOMES
GENOMES=( 1BU  ${REFBASE}/1BU_Barisic/1BUreseq_contigs_1K_cdhit.fasta
		   C3 ${REFBASE}/Acastellanii_genome_annotation/C3/Acanthamoeba_castellanii_C3.scaffolds.fa
		   Neff ${REFBASE}/Acastellanii_Neff/Acanthamoeba_castellanii_Neff.scaffolds.fa
		   ADB_Neff ${REFBASE}/AmoebaDB-67_AcastellaniiNeff/AmoebaDB-67_AcastellaniiNeff_Genome.fasta )

#PrepareReference_1st_pass
for i in ${!GENOMES[@]};
do
	echo "at $(date) call createSTARidx_noanno $i ${GENOMES[$i]}" >> $LOGFILE  
	createSTARidx_noanno "$i" ${GENOMES["$i"]}
done


for BASENAME in ${BASENAMES[@]}; 
do
	#Trimm Reads
	R1=${RAWDATA}/${BASENAME}_S[0-9][0-9]_[R]1_*.fastq.gz
	R2=${RAWDATA}/${BASENAME}_S[0-9][0-9]_[R]2_*.fastq.gz
	TrimmReads $R1 $R2 $BASENAME

	# Remove host reads
	T1=${WORKDIR}/TrimmedReads/${BASENAME}_trimmed_1.fastq.gz
	T2=${WORKDIR}/TrimmedReads/${BASENAME}_trimmed_2.fastq.gz

	filterHuman $T1 $T2 $BASENAME
	for gen_bn in ${!GENOMES[@]};
	do
		echo "at $(date) map $BASENAME to $gen_bn" >> $LOGFILE  
		mapSTAR_first $BASENAME $gen_bn
done

done

