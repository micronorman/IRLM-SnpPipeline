#!/usr/bin/env sh

# Path to Samtools
SAMTOOLS="$( exec_path samtools )"
[ $? -eq 0 ] || croak "samtools not found in PATH"

SAMTOOLS_HOME="$( dirname "$SAMTOOLS" )"
SAMTOOLS_VERSION="$("$SAMTOOLS" 2>&1 | grep ^Version | cut -f2 -d' ' )"

SAM_EXTS=('sam' 'bam' 'cram')
SAM_REGEX="\.($( join $'|' "${SAM_EXTS[@]}" )"

# Path to Bcftools
BCFTOOLS="$( exec_path bcftools )"
[ $? -eq 0 ] || croak "bcftools not found in PATH"

BCFTOOLS_VERSION="$("$BCFTOOLS" 2>&1 | grep ^Version | cut -f2 -d' ' )"

GREP_CMD="$( exec_path grep )"

# Options
SAMTOOLS_COMPRESSION_LEVEL=6
SAMTOOLS_OUTPUT_FORMAT='bam'
SAMTOOLS_PILEUP_MAX_DEPTH=8000
SAMTOOLS_SORT_MAXMEM="$(( (8 * SYSMEM) / (10 * NUM_THREADS) ))K"
SAMTOOLS_NOFILTER=false
SAMTOOLS_NOSORT=false
SAMTOOLS_EXCLUDE_BITS='0x400'

is_sam () {
	debug "$*"

	for fileName in "$@"; do
		for ext in "${SAM_EXTS[@]}"; do
			if [[ "$fileName" =~ $ext$ ]]; then
				echo "$ext"
				return 0
			fi
		done
	done

	return 1
}

find_sam () {
	debug "$*"
	local inputPrefix="$1"
	
	find_files "$inputPrefix" "$SAM_REGEX" 'FASTQ-'

	return $?
}

make_fai_index () {
	debug "$*"
	local inputPrefix="$1"
	local inputFastaFile=

	if is_fasta "$inputPrefix" > /dev/null; then
		inputFastaFile="$inputPrefix"
	else 
		inputFastaFile="$( find_fasta "$inputPrefix" )"
	fi

	[ -n "$inputFastaFile" ] || croak "Could not index sequence file for prefix $inputPrefix"

	info "Indexing fasta-file $inputFastaFile"

	exec_command $SAMTOOLS faidx "${inputFastaFile}"

	return $? 
}

make_bai_index () {
	debug "$*"
	local inputBamFile="$1"

	check_exists "$inputBamFile" "BAM"

	info "Indexing BAM-file $inputBamFile"

	exec_command $SAMTOOLS "index" "$inputBamFile"		
}

fix_mates () {
	debug "$*"
	local inputBamFile="$1"
	local outputBamFile="$2"

	check_exists "$inputBamFile" "SAM/BAM"

	exec_command $SAMTOOLS "fixmate" "$inputBamFile" "$outputBamFile"
}

#generate_pileup () {
#	debug "$*"
#	local inputSamFile="$1"
#	local outputBcfFile="$2"
#	local referenceFastaFile="$3"
#
#	check_exists "$inputSamFile" "Input alignment SAM/BAM"
#
#	
#}

call_variants_quick () {
	local inputBamFile="$1"
	local referenceFastaFile="$2"
	local outputPrefix="$3"

	#if ! is_fasta "$referenceFastaFile"; then
	#	local referenceFolder="$(dirname "$referenceFastaFile")"
	#	referenceFastaFile="$(find "$referenceFolder/" | grep -E "$referenceFastaFile\.(fa|fna|fasta)$" | head -n 1 )"
	#fi

	check_exists "$referenceFastaFile" 'reference FASTA'

	samtoolsPileupArgs=(mpileup -uf "$referenceFastaFile" "$inputBamFile")

	# Versions 0.1.18+ uses legacy bcftools viewer to call variants. Version >1.0 
	# uses the separate bcftools call program instead
	bcftoolsArgs=(call -mvo "$outputPrefix.vcf")

	finalExec=($SAMTOOLS "${samtoolsPileupArgs[@]}" "|" \
	           $BCFTOOLS "${bcftoolsArgs[@]}" \
	)

	exec_command "${finalExec[@]}"

	return $?
}

filter_alignments () {
	local inputSamFil="$1"
	
}

filter_variants () {
	local inputVCFFile="$1"
	local outputFile="$2"
	local minDepth="$3"
	local minMapQ="$4"
	local minAltPerc="$5"

	check_exists "$inputVCFFile" 'VCF'

	if [ -z "$outputFile" ]; then
		exec 5>&1
	else
		exec 5>"$outputFile"
	fi

	if [ -z "$minDepth" ]; then minDepth=5; fi
	if [ -z "$minMapQ" ]; then minMapQ=45; fi
	if [ -z "$minAltPerc" ]; then minAltPerc=0.7; fi

	exec_command "$SCRIPT_HOME/filter_snps.pl" -d "$minDepth" -m "$minMapQ" -p "$minAltPerc" "$inputVCFFile" 1>&5

	exec 5>&-
}

generate_mapping_stats () {
	local inputBamFile="$1"
	local referencePrefix="$2"
	local outputFile="$3"
	local referenceFastaFile=

	if is_fasta "$referencePrefix" 1> /dev/null; then
		referenceFastaFile="$referencePrefix"
	else
		referenceFastaFile="$( find_fasta "$referencePrefix" )"	
	fi

	check_exists "$inputBamFile" 'input BAM/SAM'
	#check_exists "$referenceFastaFile" 'reference'

	if [ -z "$outputFile" ]; then
		exec 5>&1
	else
		exec 5>"$outputFile"
	fi

	#exec_command $SAMTOOLS "stats" "-r $referenceFastaFile" -F "$SAMTOOLS_EXCLUDE_BITS" "$inputBamFile" 1>&5
	
	exec_command $SAMTOOLS "stats" -m '0.98' -F "$SAMTOOLS_EXCLUDE_BITS" "$inputBamFile" 1>&5

	local exitStatus=$?

	exec 5>&-

	return $exitStatus
}

stats_plot_bamstats () {
	local inputStatsFile="$1"
	local outputPrefix="$2"

	check_exists "$inputStatsFile" 'Samtools stats'

	if [ -z "$outputPrefix" ]; then
		outputPrefix="$( basename "$inputStatsFile" ).plots"
	fi

	exec_command "$SAMTOOLS_HOME/plot-bamstats" "-p $outputPrefix" "$inputStatsFile"
}

stats_coverage () {
	local inputStatsFile="$1"
	local outputFile="$2"

	check_exists "$inputStatsFile" 'Samtools stats'

	if [ -z "$outputFile" ]; then
		exec 5>&1
	else
		exec 5>"$outputFile"
	fi

	exec_command $GREP_CMD "^COV" "$inputStatsFile" "|" cut "-f" "2-" 1>&5

	exec 5>&-
}


calculate_coverage () {
	local inputBamFile="$1"
	local outputFile="$2"

	check_exists "$inputBamFile" 'input BAM/SAM'

	if [ -z "$outputFile" ]; then
		exec 5>&1
	else
		exec 5>"$outputFile"
	fi

	exec_command "$SCRIPT_HOME/calc_avg_bam_cov.pl" "$inputBamFile" 1>&5  

	exec 5>&-
}

sam2sortedbam_quick() {
	local inputSamFile="$1"
	local outputPrefix="$2"

	check_exists "$inputSamFile" 'SAM'

	if [ -z "$outputPrefix" ]; then
		outputPrefix="$( basename "$inputSamFile" .sam )"
	fi

	samtoolsSortArgs=( -@ "$NUM_THREADS" \
	                   -m "$SAMTOOLS_SORT_MAXMEM" \
	                   -l 0 \
	                   -O "bam" \
	                   -T "_tmp.$( basename "$inputSamFile" .sam)" \
	                   -o "${outputPrefix}.bam" \
	)

	exec_command $SAMTOOLS "sort" "${samtoolsSortArgs[@]}" "$inputSamFile" 
}

coordinatesort_sam() {
	local inputPrefix="$1"
	local outputPrefix="$2"

	local inputSamFile=
	local samExt="$( is_sam "$inputPrefix" )"
	
	if [ $? -eq 0 ];then
		inputSamFile="$inputPrefix";
		inputPrefix="$( basename "$inputPrefix" ".$samExt" )"
	else
		inputSamFile="$( find_sam "$inputPrefix" )"
	fi 

	if [ -z "$inputSamFile" ]; then
		croak "No input SAM file found for prefix $inputPrefix"
	fi

	if [ -z "$outputPrefix" ]; then
		outputPrefix="$inputPrefix"
	fi

	local outputSamFile="$outputPrefix.$SAMTOOLS_OUTPUT_FORMAT"

	# NOTE:Only works with samtools-1.0+ at this point	
	local samtoolsArgs=( sort \
	                     -@ "$NUM_THREADS" \
	                     -m "$SAMTOOLS_SORT_MAXMEM" \
	                     -l "$SAMTOOLS_COMPRESSION_LEVEL" \
	                     -O "$SAMTOOLS_OUTPUT_FORMAT"
	                     -T "$outputPrefix" \
	)

	exec_command $SAMTOOLS "${samtoolsArgs[@]}" "$inputSamFile" ">$outputSamFile"
}

store_unmapped_reads() {
	local inputPrefix="$1"
	local outputPrefix="$2"

	local inputSamFile=
	local samExt="$( is_sam "$inputPrefix" )"

	if [ $? -eq 0 ];then
		inputSamFile="$inputPrefix";
		inputPrefix="$( basename "$inputPrefix" ".$samExt" )"
	else
		inputSamFile="$( find_sam "$inputPrefix" )"
	fi 

	if [ -z "$inputSamFile" ]; then
		croak "No input SAM file found for prefix $inputPrefix"
	fi

	if [ -z "$outputPrefix" ]; then
		outputPrefix="$inputPrefix"
	fi

	local samtoolsFastqArgs=( fastq \
							  -f '0xC' \
							  -c "$GZIP_COMPRESS_LEVEL" \
							  -1 "$outputPrefix.unmapped_R1.fastq.gz" \
							  -2 "$outputPrefix.unmapped_R2.fastq.gz" \
	)

	exec_command $SAMTOOLS "${samtoolsFastqArgs[@]}" "$inputSamFile"
}

finalize_sam() {
	local inputPrefix="$1"
	local outputPrefix="$2"
	local referencePrefix="$3"
	local finalOutputFormat="$4"

	local referenceFastaFile=

	if is_fasta "$referencePrefix" 1> /dev/null; then
		referenceFastaFile="$referencePrefix"
	else
		referenceFastaFile="$( find_fasta "$referencePrefix" )"	
	fi

	if [ -z "$finalOutputFormat" ]; then
		finalOutputFormat="$SAMTOOLS_OUTPUT_FORMAT"
	fi

	local inputSamFile=
	local samExt="$( is_sam "$inputPrefix" )"
	
	if [ $? -eq 0 ];then
		inputSamFile="$inputPrefix";
		inputPrefix="$( basename "$inputPrefix" ".$samExt" )"
	else
		inputSamFile="$( find_sam "$inputPrefix" )"
	fi 

	if [ -z "$inputSamFile" ]; then
		croak "No input SAM file found for prefix $inputPrefix"
	fi

	if [ -z "$outputPrefix" ]; then
		outputPrefix="$inputPrefix"
	fi

	local outputSamFile="$outputPrefix.$finalOutputFormat"
	local samtoolsViewSwitch='-u'

	if [ "$finalOutputFormat" = 'cram' ]; then
		samtoolsViewSwitch='-C'
	fi

	local samtoolsFilterArgs=( view "${samtoolsViewSwitch}" \
                               -F "$SAMTOOLS_EXCLUDE_BITS" \
	                           -q "$READS_MINMAPQ" \
	                           -T "$referenceFastaFile" \
	)

	local samtoolsSortArgs=( sort \
	                         -@ "$NUM_THREADS" \
	                         -m "$SAMTOOLS_SORT_MAXMEM" \
	                         -l "$SAMTOOLS_COMPRESSION_LEVEL" \
	                         -O "$finalOutputFormat"
	                         -T "$outputPrefix" \
	)

	exec_command $SAMTOOLS "${samtoolsFilterArgs[@]}" "$inputSamFile" "|" $SAMTOOLS "${samtoolsSortArgs[@]}" ">$outputSamFile"
}

#namesort_bam () {
#	local inputSamFile="$1"
#}


sam2unpairedfq() {
	local inputSamFile="$1"
	local outputPrefix="$2"
	
	check_exists "$inputSamFile" 'input SAM'

	if [ -z "$outputPrefix" ]; then
		outputPrefix="$(basename "$inputSamFile" .sam)"
	fi

	init_work_folder "$inputSamFile"

	inputSamFile="$(basename "$inputSamFile")"
	local outputFolder="$(dirname "$outputPrefix")"
	outputPrefix="out/$(basename "$outputPrefix")"
	
	local samtoolsViewArgs=(view -@ "$NUM_THREADS" -Su -F 0x2  "$inputSamFile")
	local samtoolsBam2FqArgs=(bam2fq - ">" "${outputPrefix}_unpaired.fq")

	local finalArgs=( $SAMTOOLS "${samtoolsViewArgs[@]}" "|"\
	                  $SAMTOOLS "${samtoolsBam2FqArgs[@]}"\
	)

	exec_command "${finalArgs[@]}"

	close_work_folder "$outputFolder"
}

debug "Plugin $(basename "${BASH_SOURCE[0]}") loaded"
debug "Samtools version: $SAMTOOLS_VERSION"
debug "Bcftools version: $BCFTOOLS_VERSION"

