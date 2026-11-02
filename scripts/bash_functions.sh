function createSTARidx_noanno () {
	# args: genome name and path to genome fasta (can end with fna, fa, fasta, needs to be unzipped)
	echo "$FUNCNAME on $1 $2 at "`date` >> $LOGFILE
    local refgen="$2"
    local gen_bn="$1"
	#local gen_bn=$( basename "${refgen%\.f[na][.ast]*}" )
    #local gen_dir=$(dirname "$refgen")
	if [[ $refgen =~ ".gz$" ]]; then
   		echo "The genome fastas need to be unzipped for STAR index generation!" >&2
		exit 1
	fi
    if [ -d ${WORKDIR}/Indeces/${gen_bn}_STAR1 ]; 
    then
        echo "${WORKDIR}/Indeces/${gen_bn}_STAR1 already exists, do nothing" >> $LOGFILE
        return
    fi
    mkdir -p ${WORKDIR}/Indeces/${gen_bn}_STAR1
    echo  "STAR --runThreadN 5 --runMode genomeGenerate --genomeSAindexNbases 11 --genomeDir ${WORKDIR}/Indeces/${gen_bn}_STAR1 --genomeFastaFiles $refgen" >> $LOGFILE
    STAR --runThreadN 5 --runMode genomeGenerate --genomeSAindexNbases 11 \
    --genomeDir ${WORKDIR}/Indeces/${gen_bn}_STAR1 --genomeFastaFiles $refgen 2> ${WORKDIR}/LOGS/${gen_bn}_star_idx.err.log
    return 
}

function TrimmReads () {
	# args: read1 read2 basename
	echo "$FUNCNAME on $1 $2 at "`date` >> $LOGFILE
	mkdir -p ${WORKDIR}/TrimmedReads
	local BN="$3"
	if [[ -f ${WORKDIR}/TrimmedReads/${BN}_trimmed_1.fastq.gz && -f ${WORKDIR}/TrimmedReads/${BN}_trimmed_2.fastq.gz ]]; 
	then
		echo "trimmed fastq already exist" >> $LOGFILE
		return
	fi
	fastp -l 40 -y -Y 30 -M 15 --cut_front --cut_tail --in1 $1 --in2 $2 \
	--out1 ${WORKDIR}/TrimmedReads/${BN}_trimmed_1.fastq.gz --out2 ${WORKDIR}/TrimmedReads/${BN}_trimmed_2.fastq.gz \
	--html ${WORKDIR}/Reports/${BN}_fastp.html --json ${WORKDIR}/LOGS/${BN}_fastp.json \
	2>> ${WORKDIR}/Reports/${BN}_fastp.err.log >> ${WORKDIR}/LOGS/${BN}_fastp.log 
	ES=$?
	echo finished at `date` with exit state $ES >> $LOGFILE
	return
}

function filterHuman () {
	# args: read1 read2 basename
	echo "$FUNCNAME on $1 $2 with BN $3 at "`date` >> $LOGFILE
	local BN=$3
	#READ1=${WORKDIR}/TrimmedReads/${BN}_trimmed_1.fastq.gz
	#READ2=${WORKDIR}/TrimmedReads/${BN}_trimmed_2.fastq.gz
	mkdir -p ${WORKDIR}/FilteredReads
	OUT=${WORKDIR}/FilteredReads/${BN}_filtered
	if [[ -f ${WORKDIR}/FilteredReads/${BN}_filtered_1.fastq.gz && -f ${WORKDIR}/FilteredReads/${BN}_filtered_2.fastq.gz ]]; 
	then
		echo "filtered fastq already exist" >> $LOGFILE
		return
	fi

	bowtie2 -p 5 -x $BT2IDX -1 $1 -2 $2 --local --un-conc-gz ${OUT}  > /dev/null 2> ${OUT}_bt2.err.log
	ES=$?
	echo finished at `date` with exit state $ES >> $LOGFILE
	mv ${OUT}.1 ${OUT}_1.fastq.gz
	mv ${OUT}.2 ${OUT}_2.fastq.gz
}


function mapSTAR_first () {
	# read bn $1 genome bn $2 
    echo "FUNCTION: $FUNCNAME on $1 $2 at "`date` >> $LOGFILE 
    mkdir -p ${WORKDIR}/STARMap1
	local ERRLOG=${WORKDIR}/LOGS/star_err.log
	local BN=$1
	local GEN_BN=$2
	#local GEN_BNS=${@}
	# get rid of first arg (BN)
	#unset GEN_BNS[0]
    #for GEN_BN in "${GEN_BNS[@]}"
    #do
	if [ -f ${WORKDIR}/STARMap1/${BN}_${GEN_BN}"Aligned.sortedByCoord.out.bam" ]; 
	then
		echo "sorted bam already exist" >> $LOGFILE
		return
	fi
	local GEN_DIR=${WORKDIR}/Indeces/${GEN_BN}_STAR1
	local R1=${WORKDIR}/FilteredReads/${BN}_filtered_1.fastq.gz
	local R2=${WORKDIR}/FilteredReads/${BN}_filtered_2.fastq.gz
	echo "mapping: $1 to $GEN_BN at "`date` >> $LOGFILE 	
	echo STAR --runMode alignReads \
								--genomeDir $GEN_DIR \
								--runThreadN 10 \
								--readFilesCommand zcat \
								--readFilesIn $R1 $R2 \
								--outFileNamePrefix ${WORKDIR}/STARMap1/${BN}_${GEN_BN} \
								--outSAMtype BAM Unsorted >> $LOGFILE
								#--outSAMmapqUnique Integer0to255 \
								#--quantMode TranscriptomeSAM >> $LOGFILE
	STAR --runMode alignReads \
							--genomeDir $GEN_DIR \
							--runThreadN 10 \
							--readFilesCommand zcat \
							--readFilesIn $R1 $R2 \
							--outFileNamePrefix ${WORKDIR}/STARMap1/${BN}_${GEN_BN} \
							--outSAMtype BAM Unsorted 2>> $ERRLOG >> $LOGFILE
							#--outSAMmapqUnique Integer0to255 \
							#--quantMode TranscriptomeSAM >> $LOGFILE
	ES=$?
	echo finished mapping at `date` with exit state $ES >> $LOGFILE
	echo "start samtools sort ${BN}_${GEN_BN} at $(date)"	>> $LOGFILE
	echo samtools sort --threads 4 -m 5G ${WORKDIR}/STARMap1/${BN}_${GEN_BN}Aligned.out.bam -o ${WORKDIR}/STARMap1/${BN}_${GEN_BN}Aligned.sortedByCoord.out.bam >> $LOGFILE
	#attention no .bam for output file -m increases memory to reduce the number of intermediate files produced by samtools that throw an error message
	samtools sort --threads 4 -m 5G ${WORKDIR}/STARMap1/${BN}_${GEN_BN}Aligned.out.bam -o ${WORKDIR}/STARMap1/${BN}_${GEN_BN}Aligned.sortedByCoord.out.bam 2>> $ERRLOG
    rm -f ${WORKDIR}/STARMap1/${BN}_${GEN_BN}Aligned.out.bam
	samtools index ${WORKDIR}/STARMap1/${BN}_${GEN_BN}Aligned.sortedByCoord.out.bam
	#done
    #on the fly sorting needs far too much RAM 		
    #doesnt work with compressed fastq files also not with the corresponding option
    #do NOT USE the sorting option. needs too much memory! --outSAMtype BAM SortedByCoordinate output sorted by coordinate
}

function mapSTAR_main () {
	# read bn $1 genome index $2 
    echo "FUNCTION: $FUNCNAME on $1 with $2 at "`date` >> $LOGFILE 
    mkdir -p ${WORKDIR}/STARMap1
	local ERRLOG=${WORKDIR}/LOGS/star_err.log
	local BN=$1
	local GEN_IDX=$2
	#local GEN_BNS=${@}
	# get rid of first arg (BN)
	#unset GEN_BNS[0]
    #for GEN_BN in "${GEN_BNS[@]}"
    #do
	if [ -f ${WORKDIR}/STARMap_main/${BN}_"Aligned.sortedByCoord.out.bam" ]; 
	then
		echo "sorted bam already exist" >> $LOGFILE
		return
	fi
	local GEN_DIR=${WORKDIR}/Indeces/$2
	local R1=${WORKDIR}/FilteredReads/${BN}_filtered_1.fastq.gz
	local R2=${WORKDIR}/FilteredReads/${BN}_filtered_2.fastq.gz
	echo "mapping: $1 to $GEN_DIR at "`date` >> $LOGFILE 	
	echo STAR --runMode alignReads \
								--genomeDir $GEN_DIR \
								--runThreadN 10 \
								--readFilesCommand zcat \
								--readFilesIn $R1 $R2 \
								--outFileNamePrefix ${WORKDIR}/STARMap_main/${BN}_ \
								--outSAMtype BAM Unsorted >> $LOGFILE
								#--outSAMmapqUnique Integer0to255 \
								#--quantMode TranscriptomeSAM >> $LOGFILE
	STAR --runMode alignReads \
							--genomeDir $GEN_DIR \
							--runThreadN 10 \
							--readFilesCommand zcat \
							--readFilesIn $R1 $R2 \
							--outFileNamePrefix ${WORKDIR}/STARMap_main/${BN}_ \
							--outSAMtype BAM Unsorted 2>> $ERRLOG >> $LOGFILE
							#--outSAMmapqUnique Integer0to255 \
							#--quantMode TranscriptomeSAM >> $LOGFILE
	ES=$?
	echo finished mapping at `date` with exit state $ES >> $LOGFILE
	echo "start samtools sort ${BN}_Aligned.out.bam at $(date)"	>> $LOGFILE
	echo samtools sort --threads 4 -m 5G ${WORKDIR}/STARMap_main/${BN}_Aligned.out.bam -o ${WORKDIR}/STARMap_main/${BN}_Aligned.sortedByCoord.out.bam >> $LOGFILE
	#attention no .bam for output file -m increases memory to reduce the number of intermediate files produced by samtools that throw an error message
	samtools sort --threads 4 -m 5G ${WORKDIR}/STARMap_main/${BN}_Aligned.out.bam -o ${WORKDIR}/STARMap_main/${BN}_Aligned.sortedByCoord.out.bam 2>> $ERRLOG
    #rm -f ${WORKDIR}/STARMap_main/${BN}_Aligned.out.bam
	samtools index ${WORKDIR}/STARMap_main/${BN}_Aligned.sortedByCoord.out.bam
	#done
    #on the fly sorting needs far too much RAM 		
    #doesnt work with compressed fastq files also not with the corresponding option
    #do NOT USE the sorting option. needs too much memory! --outSAMtype BAM SortedByCoordinate output sorted by coordinate
}


function RNAseqPicardQC () {
	# gets basename and mapped RNA file
	# create ref_flat annotation file from gff3 using gff3ToGenePred (mamba install ucsc-gff3togenepred)
	# gff3ToGenePred  -bad=Acanthamoeba_castellanii_C3_genpred.bad Acanthamoeba_castellanii_C3.gff3 Acanthamoeba_castellanii_C3.genepred
	# awk -v FS="\t" -v OFS="\t" '{print $12,$1, $2, $3, $4, $5, $6, $7, $8 ,$9, $10}' ${REFBASE}/Acastellanii_C3/Acanthamoeba_castellanii_C3.genepred > ${REFBASE}/Acastellanii_C3/Acanthamoeba_castellanii_C3.ref_flat
	# for ribosomal interval list create genome dictionary and rRNA locations:
	# samtools dict Acanthamoeba_castellanii_C3.scaffolds.fa > Acanthamoeba_castellanii_C3.scaffolds.dict
	# cat Acanthamoeba_castellanii_C3.scaffolds.dict <(awk -v FS="\t" -v OFS="\t" '$3 == "rRNA" {print $1,$5,$5,$7,$9}' companion_C3/scaffold.out.gff3 ) > Acanthamoeba_castellanii_C3_rRNA.interval_list
	echo "$FUNCNAME on file $1 at "`date` >> $LOGFILE
	local FILEBASE=$(basename $1 .bam)
	FILEBASE=${FILEBASE%%_Aligned.*}
	java -Xmx6g -jar $PICARD CollectRnaSeqMetrics \
		I=$1 O=${WORKDIR}/LOGS/${FILEBASE}_rnaseq_picard.out \
		CHART=${WORKDIR}/Reports/${FILEBASE}_rnaseq_picard.pdf \
		REF_FLAT=${GENOME_REF_FLAT} \
		STRAND=SECOND_READ_TRANSCRIPTION_STRAND \
		RIBOSOMAL_INTERVALS=${GENOME_RIBOS} \
		VALIDATION_STRINGENCY=LENIENT VERBOSITY=WARNING > ${WORKDIR}/LOGS/${FILEBASE}_RNAseq_picard.log 2>&1
}

function Count_featureCounts () {
	echo "$FUNCNAME on file $1 for $2 at "`date` >> $LOGFILE
	local BN=$2
	mkdir -p ${WORKDIR}/FeatureCounts
	mkdir -p "${WORKDIR}/${BN}_tmp"
	featureCounts -T 5 -d 30 -s 2 -Q 1 -C -p --countReadPairs -B \
				-g gene_id -F GTF -a $GENOME_GTF \
				--minOverlap 10 --largestOverlap \
				-o ${WORKDIR}/FeatureCounts/${BN}_fc.txt \
				--tmpDir ${WORKDIR}/${BN}_tmp $1 2> ${WORKDIR}/LOGS/${BN}_fc.log
	echo "finished $FUNCNAME on file $1 for $2 at "`date` >> $LOGFILE
	rm -rf "${WORKDIR}/${BN}_tmp"
}

