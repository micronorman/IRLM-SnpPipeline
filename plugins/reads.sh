#!/usr/bin/env sh

# Allowed file extensions for FastA and FastQ files
FASTQ_WHITELIST=('fq' 'fastq')
FASTA_WHITELIST=('fa' 'fna' 'fas' 'fasta')

FASTQ_REGEX="\.($( join $'|' "${FASTQ_WHITELIST[@]}" ))(\.gz)?"
FASTA_REGEX="\.($( join $'|' "${FASTA_WHITELIST[@]}" ))(\.gz)?"

# The default minimum mapping- and base- qualities. Only applies to the fastq-format
if [ -z "$READS_MINMAPQ" ]; then READS_MINMAPQ=0; fi
if [ -z "$READS_MINBAQ" ];then READS_MINBAQ=13; fi

unpack_pe_fastq_reads () {
	debug "$*"
	local inputPrefix="$1"
	local outputPrefix="$2"

	local inputFiles=( $( find_fastq "$inputPrefix" ) )

	init_work_folder "${inputFiles[@]}"

	# Strip directories from input file names as they are now linked in the working folder
	inputFiles=( $( strip_paths "${inputFiles[@]}" ) )
	numInputFiles="${#inputFiles[@]}"

	if [ "$numInputFiles" -ne 2 ]; then
		croak "Wrong number of input files found ($numInputFiles) for prefix $inputPrefix"	
	fi

	local outputFolder="$(dirname "$outputPrefix")"
	outputPrefix="$(basename "$outputPrefix")"

	# Unzip reads for sickle if they are zipped
	if is_gzipped "${inputFiles[0]}"; then
		gzip_unpack ${inputFiles[0]} "out/${outputPrefix}_1.fq"
		gzip_unpack ${inputFiles[1]} "out/${outputPrefix}_2.fq"
	else
		ln -s ${inputFiles[0]} "out/${outputPrefix}_1.fq"
		ln -s ${inputFiles[1]} "out/${outputPrefix}_2.fq"
	fi

	close_work_folder "$outputFolder"

	echo "${outputPrefix}_1.fq"
	echo "${outputPrefix}_2.fq"
 
	debug "End"
}

find_fasta () {
	debug "$*"
	local inputPrefix="$1"
	
	find_files "$inputPrefix" "$FASTA_REGEX" 'FASTA-'

	return $?
}

find_fastq () {
	debug "$*"
	local inputPrefix="$1"
	
	find_files "$inputPrefix" "$FASTQ_REGEX" 'FASTQ-'

	return $?
}

is_fastq () {
	debug "$*"
	for fileName in "$@"; do
		for ext in "${FASTQ_WHITELIST[@]}"; do
			if [[ "$fileName" =~ $ext$ ]]; then
				echo "$ext"
				return 0
			fi
		done
	done

	return 1
}

is_fasta () {
	for fileName in "$@"; do
		if [ ! -e "$fileName" ]; then return 1; fi

		for ext in "${FASTA_WHITELIST[@]}"; do
			if [[ "$fileName" =~ $ext$ ]]; then
				echo "$ext"
				return 0
			fi
		done
	done

	return 1
}

# Report loading of module to the debugger
debug "$(basename "${BASH_SOURCE[0]}") plugin loaded"
debug "FASTA_REGEX: $FASTA_REGEX"
debug "FASTQ_REGEX: $FASTQ_REGEX"
