#!/bin/bash
#----------
# author: Lukas Endler
# authored: Monday, June 3, 2024 @ 14:19:00 PM
# Time-stamp: <Friday, May 24, 2024 @ 12:33:06 PM lukasendler>
#--------------

WORKDIR=$(pwd)
echo WORKDIR=$WORKDIR
SCRIPTDIR=$(dirname -- "$0")
echo SCRIPTDIR=$SCRIPTDIR

BASEDIR=/shared_data/shared_data_BACKUP/experimental_data/P_2024_ACRS_AcanthamoebaRNASeq_LE
RAWDATA=${BASEDIR}/Reads/
REFBASE=${BASEDIR}/References/Acastellanii_genome_annotation/C3/
BT2IDX="/shared_data/shared_data_BACKUP/auxiliary_data/indices/Human/GRCh38.p14/GRCh38_noalt_as_bowtie2/GRCh38_noalt_as"
SCRIPTS=/${BASEDIR}/scripts/
PICARD="${HOME}/Tools/picard.jar"
LOGFILE=${WORKDIR}/"logfile.log"
ERRLOG=${WORKDIR}/"logfile.err.log"

source "${SCRIPTDIR}/bash_functions.sh"

mkdir -p ${WORKDIR}/LOGS
mkdir -p ${WORKDIR}/Reports
BASENAMES=$@

REFGENOME=${REFBASE}/Acanthamoeba_castellanii_C3.scaffolds.fa
GENOME_GFF=${REFBASE}/Acanthamoeba_castellanii_C3.gff3
GENOME_GTF=${REFBASE}/Acanthamoeba_castellanii_C3_gffread.gtf
GENOME_REF_FLAT=${REFBASE}/Acanthamoeba_castellanii_C3.ref_flat
GENOME_RIBOS=${REFBASE}/Acanthamoeba_castellanii_C3_rRNA.interval_list
echo "at $(date) call createSTARidx_noanno C3 $REFGENOME $REF_GFF" >> $LOGFILE  
createSTARidx_anno C3 $REFGENOME $REF_GFF

for BASENAME in ${BASENAMES[@]}; 
do
        echo "at $(date) map $BASENAME to C3" >> $LOGFILE  
        mapSTAR_main $BASENAME "C3_STARanno"
        echo "at $(date) featurecounts $BASENAME to C3" >> $LOGFILE 
        RNAseqPicardQC ${WORKDIR}/STARMap_main/${BASENAME}_Aligned.sortedByCoord.out.bam  
        Count_featureCounts  ${WORKDIR}/STARMap_main/${BASENAME}_Aligned.sortedByCoord.out.bam $BASENAME 
done


