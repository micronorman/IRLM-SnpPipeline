#!/usr/bin/env bash

PICARD="$( exec_path picard )"

[ $? -eq 0 ] || croak "Picard tools executable not found in PATH. Please install picard tools (https://broadinstitute.github.io/picard) to use this module"

PICARD_VALIDATION_STRINGENCY='STRICT'

# NOTE: Compression level set to zero, since most Picard operations are performed on temporary files
PICARD_COMPRESSION_LEVEL=0

clean_sam () {
	local inputFile="$1"
	local outputFile="$2"

	shift 2

	check_exists "$inputFile" "SAM/BAM"

	if [ -z "$outputFile" ];then
		outputFile="$( relabel_filename "$inputFile" "Cleaned" )"
	fi

	local cmdArgs=( INPUT="$inputFile" \
	                OUTPUT="$outputFile" \
	)

	exec_picard "CleanSam" "${cmdArgs[@]}" "$@"
}

mark_duplicates () {
	local inputFile="$1"
	local outputFile="$2"
	
	shift 2

	check_exists "$inputFile" "SAM/BAM"

	if [ -z "$outputFile" ]; then
		outputFile="$( relabel_filename "$inputFile" "MarkedDuplicates" )"
	fi

	local cmdArgs=( INPUT="$inputFile" \
	                OUTPUT="$outputFile" \
	                METRICS_FILE="picard.MarkDuplicates.metrics" \
	)

	exec_picard "MarkDuplicates" "${cmdArgs[@]}" "$@"
}

make_dict_file () {
	inputPrefix="$1"
	shift

	local inputFastaFile=

	if is_fasta "$inputPrefix"; then
		inputFastaFile="$inputPrefix"
		inputPrefix="$( basename "$inputPrefix" ."$( is_fasta "$inputPrefix" )" )"
	else 
		inputFastaFile="$( find_fasta "$inputPrefix" )"
	fi

	local picardArgs=( REFERENCE="$inputFastaFile" \
	                   OUTPUT="$inputPrefix.dict" \
	)

	exec_picard "CreateSequenceDictionary" "${picardArgs[@]}" "$@"

	return $?	
}

exec_picard () {
	local picardModule="$1"
	shift

	picardArgs=( "$@" )

	if $DEBUG; then
		picardArgs=( "${picardArgs[@]}" VERBOSITY="DEBUG" ) 
	fi

	if [ -n "$CURR_WORK_DIR" ]; then
		picardArgs=( "${picardArgs[@]}" TMP_DIR="$CURR_WORK_DIR" )
	fi

	local picardArgs=( "${picardArgs[@]}" \
	             VALIDATION_STRINGENCY="$PICARD_VALIDATION_STRINGENCY" \
	             COMPRESSION_LEVEL="$PICARD_COMPRESSION_LEVEL" \
	)

	exec_command "$PICARD" "$picardModule" "${picardArgs[@]}"
	
	return $?
}

debug "Plugin $(basename "${BASH_SOURCE[0]}") loaded"
