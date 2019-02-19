#!/usr/bin/env bash

# File locations
TRIMMOMATIC="$( exec_path trimmomatic )"

if [ -z "$TRIMMOMATIC" ]; then
	croak "Trimmomatic not found. Please install Trimmomatic (http://www.usadellab.org/cms/index.php?page=trimmomatic) to use this module"
fi

TRIMMOMATIC_HOME="$(dirname $TRIMMOMATIC)/.."
TRIMMOMATIC_ADAPTERS_HOME="$TRIMMOMATIC_HOME/share/trimmomatic/adapters"
TRIMMOMATIC_VERSION="$($TRIMMOMATIC -version)"
# TRIMMOMATIC_JAR="$( $REALPATH "$TRIMMOMATIC_HOME/trimmomatic.jar" )"

### General Trimming options ###

# Apadter clipping
TRIMMOMATIC_CLIPADAPT=true
TRIMMOMATIC_CLIPADAPT_SEEDMM=2
TRIMMOMATIC_CLIPADAPT_PALINTHR=24
TRIMMOMATIC_CLIPADAPT_SIMPLETHR=8

# Filter short and low-quality reads
TRIMMOMATIC_FILTER=true
TRIMMOMATIC_FILTER_MINLEN=36
TRIMMOMATIC_FILTER_AVGQUAL=13

# Keep unpaired reads
TRIMMOMATIC_KEEP_UNPAIRED=false

# Keep trim log
TRIMMOMATIC_KEEP_TRIMLOG=false

# Trimming options for pre-processing
TRIMMOMATIC_ENDTRIM=true
TRIMMOMATIC_ENDTRIM_LEADQUAL=3
TRIMMOMATIC_ENDTRIM_TRAILQUAL=3
TRIMMOMATIC_TOPHRED33=true

# Trimming options for read trimming
TRIMMOMATIC_SLIDING_WINSIZE=4
TRIMMOMATIC_SLIDING_AVGQUAL=13

pre_process_pe_reads () {
	# Clips adapters and leading and trailing low-quality bases. This function only uses end trimming
	# and does not filter on average read quality

	local inputPrefix="$1"
	local outputPrefix="$2"
	local adaptersPrefix="$3"

	inputFiles=( $( find_fastq "$inputPrefix" ) )

	if [ -z "$adaptersPrefix" ]; then
		adaptersPrefix="NexteraPE"
	fi

	local adaptersFile="$TRIMMOMATIC_ADAPTERS_HOME/${adaptersPrefix}-PE.fa"

	check_exists "$adaptersFile" 'Library adapters'

	init_work_folder "${inputFiles[@]}"

	inputFiles=( $( strip_paths "${inputFiles[@]}" ) )
	
	# Check that we have two input files
	if [ "${#inputFiles[@]}" -lt 2 ]; then
		croak "Not enough input-files found for paired-end trimming using input prefix $inputPrefix"
	fi
	
	# Define output file names
	local outputFolder="$( dirname "$outputPrefix" )"
	
	outputPrefix="$( basename "$outputPrefix" )"

	local outTrimLog="${outputPrefix}.trimlog.txt"
	local outPE1="${outputPrefix}_1.fq"
	local outPE2="${outputPrefix}_2.fq"
	local outUP1="${outputPrefix}_1up.fq"
	local outUP2="${outputPrefix}_2up.fq"

	# Build command line arguments
	local trimmomaticArgs=( "PE" -threads "$NUM_THREADS" -trimlog "$outTrimLog" )
	

	trimmomaticArgs=( "${trimmomaticArgs[@]}" \
	                  "${inputFiles[@]}" \
	                  "$outPE1" \
	                  "$outUP1" \
	                  "$outPE2" \
	                  "$outUP2" \
	)

	if $TRIMMOMATIC_CLIPADAPT; then
		local illuminaClipFields=( "$adaptersFile" \
		                           "$TRIMMOMATIC_CLIPADAPT_SEEDMM" \
		                           "$TRIMMOMATIC_CLIPADAPT_PALINTHR" \
		                           "$TRIMMOMATIC_CLIPADAPT_SIMPLETHR" \
		)

		local illuminaClipString="$( join ":" "${illuminaClipFields[@]}" )"

		trimmomaticArgs=( "${trimmomaticArgs[@]}" "ILLUMINACLIP:$illuminaClipString" )
	fi

	if $TRIMMOMATIC_ENDTRIM; then
		trimmomaticArgs=( "${trimmomaticArgs[@]}" \
	                          "LEADING:$TRIMMOMATIC_ENDTRIM_LEADQUAL" \
	                          "TRAILING:$TRIMMOMATIC_ENDTRIM_TRAILQUAL" \
		)
	fi

	if $TRIMMOMATIC_FILTER; then
		trimmomaticArgs=( "${trimmomaticArgs[@]}" \
	                  "MINLEN:$TRIMMOMATIC_FILTER_MINLEN" \
		)
	fi

	if $TRIMMOMATIC_TOPHRED33; then
		trimmomaticArgs=( "${trimmomaticArgs[@]}" "TOPHRED33" )
	fi

	exec_command "$TRIMMOMATIC" "${trimmomaticArgs[@]}"

	if [ $? -eq 0 ]; then 
		gzip_pack "$outPE1" "out/$outPE1.gz"
		gzip_pack "$outPE2" "out/$outPE2.gz"
	
		if $TRIMMOMATIC_KEEP_UNPAIRED; then
			gzip_pack "$outUP1" "out/$outUP1.gz"
			gzip_pack "$outUP2" "out/$outUP2.gz"
		fi
	fi

	if $TRIMMOMATIC_KEEP_TRIMLOG; then
		gzip_pack "$outTrimLog" "out/$outTrimLog"
	fi

	close_work_folder "$outputFolder"
}

pre_process_se_reads () {
	local inputPrefix="$1"
	local outputPrefix="$2"
	local adaptersPrefix="$3"

	inputFiles=( $( find_fastq "$inputPrefix" ) )

	if [ -z "$adaptersPrefix" ]; then
		adaptersPrefix="All"
	fi

	local adaptersFile="$TRIMMOMATIC_ADAPTERS_HOME/${adaptersPrefix}-SE.fa"

	check_exists "$adaptersFile" 'Library adapters'

	init_work_folder "${inputFiles[@]}"

	inputFiles=( $( strip_paths "${inputFiles[@]}" ) )

	# Define output file names
	local outputFolder="$( dirname "$outputPrefix" )"
	
	outputPrefix="$( basename "$outputPrefix" )"

	local outTrimLog="${outputPrefix}.trimlog.txt"
	local outSE="${outputPrefix}.fq"

	# Build command line arguments
	local trimmomaticArgs=( "SE" -threads "$NUM_THREADS" -trimlog "$outTrimLog" )

	trimmomaticArgs=( "${trimmomaticArgs[@]}" \
	                  "${inputFiles[@]}" \
	                  "$outSE" \
	)

	if $TRIMMOMATIC_CLIPADAPT; then
		local illuminaClipFields=( "$adaptersFile" \
		                           "$TRIMMOMATIC_CLIPADAPT_SEEDMM" \
		                           "$TRIMMOMATIC_CLIPADAPT_PALINTHR" \
		                           "$TRIMMOMATIC_CLIPADAPT_SIMPLETHR" \
		)

		local illuminaClipString="$( join ":" "${illuminaClipFields[@]}" )"

		trimmomaticArgs=( "${trimmomaticArgs[@]}" "ILLUMINACLIP:$illuminaClipString" )
	fi

	if $TRIMMOMATIC_ENDTRIM; then
		trimmomaticArgs=( "${trimmomaticArgs[@]}" \
	                          "LEADING:$TRIMMOMATIC_ENDTRIM_LEADQUAL" \
	                          "TRAILING:$TRIMMOMATIC_ENDTRIM_TRAILQUAL" \
		)
	fi

	if $TRIMMOMATIC_FILTER; then
		trimmomaticArgs=( "${trimmomaticArgs[@]}" \
	                  "MINLEN:$TRIMMOMATIC_FILTER_MINLEN" \
		)
	fi

	if $TRIMMOMATIC_TOPHRED33; then
		trimmomaticArgs=( "${trimmomaticArgs[@]}" "TOPHRED33" )
	fi

	exec_command "$TRIMMOMATIC" "${trimmomaticArgs[@]}"

	if [ $? -eq 0 ]; then 
		gzip_pack "$outSE" "out/$outSE.gz"
	else
		warn "Trimmomatic did not run succesfully"
		return 1
	fi

	if $TRIMMOMATIC_KEEP_TRIMLOG; then
		gzip_pack "$outTrimLog" "out/$outTrimLog"
	fi

	close_work_folder "$outputFolder"
}

trim_pe_reads () {
	local inputPrefix="$1"
	local outputPrefix="$2"

	inputFiles=( $( find_fastq "$inputPrefix" ) )

	init_work_folder "${inputFiles[@]}"

	inputFiles=( $( strip_paths "${inputFiles[@]}" ) )
	
	# Check that we have two input files
	if [ "${#inputFiles[@]}" -lt 2 ]; then
		croak "Not enough input-files found for paired-end trimming using input prefix $inputPrefix"
	fi
	
	# Define output file names
	local outputFolder="$( dirname "$outputPrefix" )"
	
	outputPrefix="$( basename "$outputPrefix" )"

	local outTrimLog="${outputPrefix}.trimlog.txt"
	local outPE1="${outputPrefix}_1.fq"
	local outPE2="${outputPrefix}_2.fq"
	local outUP1="${outputPrefix}_1up.fq"
	local outUP2="${outputPrefix}_2up.fq"


	# Build command line arguments
	local trimmomaticArgs=( "PE" -threads "$NUM_THREADS" -trimlog "$outTrimLog" )

	trimmomaticArgs=( "${trimmomaticArgs[@]}" \
	                  "${inputFiles[@]}" \
	                  "$outPE1" \
	                  "$outUP1" \
	                  "$outPE2" \
	                  "$outUP2" \
	)

	trimmomaticArgs=( "${trimmomaticArgs[@]}" \
			  "SLIDINGWINDOW:${TRIMMOMATIC_SLIDING_WINSIZE}:${TRIMMOMATIC_SLIDING_AVGQUAL}" \
	)

	if $TRIMMOMATIC_FILTER; then
		trimmomaticArgs=( "${trimmomaticArgs[@]}" \
	                  "MINLEN:$TRIMMOMATIC_FILTER_MINLEN" \
		          "AVGQUAL:$TRIMMOMATIC_FILTER_AVGQUAL"
		)
	fi

	exec_command "$TRIMMOMATIC" "${trimmomaticArgs[@]}"

	if [ $? -eq 0 ]; then 
		gzip_pack "$outPE1" "out/$outPE1.gz"
		gzip_pack "$outPE2" "out/$outPE2.gz"
	
		if $TRIMMOMATIC_KEEP_UNPAIRED; then
			gzip_pack "$outUP1" "out/$outUP1.gz"
			gzip_pack "$outUP2" "out/$outUP2.gz"
		fi
	else
		warn "Trimmomatic did not run succesfully"
		return 1

	fi

	if $TRIMMOMATIC_KEEP_TRIMLOG; then
		gzip_pack "$outTrimLog" "out/$outTrimLog"
	fi

	close_work_folder "$outputFolder"
}

debug "Plugin $(basename "${BASH_SOURCE[0]}") loaded"
debug "Trimmomatic version: $TRIMMOMATIC_VERSION"
debug "Trimmomatic adaptesrs path: $TRIMMOMATIC_ADAPTERS_HOME"
