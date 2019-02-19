#!/usr/bin/env sh

# Paths to the BWA mapping tool
BWA_HOME="$(dirname "$( exec_path bwa )")"

if [ -z "$BWA_HOME" ]; then
	croak "No bwa installation was found. Cannot load module"
fi

BWA_CMD="$BWA_HOME/bwa"

# BWA version control
BWA_VERSION="$($BWA_CMD 2>&1 | grep ^Version | cut -f 2 -d' ')"
BWA_MINOR_VERSION="$(echo $BWA_VERSION | cut -f 2 -d '.' )"
BWA_BUILD="$(echo $BWA_VERSION | cut -f 2 -d'-' )"

BWA_REQ_MINOR_VERSION=6

# Pre- bwa 0.6 had more files. So this is just in case they decide to change it up again
BWA_POST_6_INDEX_EXTS=( "bwt" "sa" "pac" "amb" "ann" )

if [ "$BWA_MINOR_VERSION" -lt "$BWA_REQ_MINOR_VERSION" ]; then
	warn "Wrong version of bwa ($BWA_MAJOR_VERSION). Must be newer than $BWA_MIN_VERSION"
	exit 1
fi

BWA_INDEX_EXTS=( "${BWA_POST_6_INDEX_EXTS[@]}" )
BWA_INDEX_NUMFILES="${#BWA_INDEX_EXTS[@]}"
BWA_INDEX_REGEX="\.($( join $'|' "${BWA_INDEX_EXTS[@]}" ))"

# BWA default parameters
BWA_QUERY_INTERLEAVED=false
BWA_MAPPING_ENGINE='mem'

build_reference_index () {
	local inputPrefix="$1"
	local outputPrefix="$2"
	local exitStatus=0

	if [[ $# -lt 1 || $# -gt 2 ]]; then
		croak "Argument error"
	fi

	local inputFastaFile=

	if is_fasta $inputPrefix; then
		inputFastaFile="$inputPrefix"
	else 
		inputFastaFile="$( find_fasta "$inputPrefix" )"
	fi

	if [ -z "$inputFastaFile" ]; then
		croak "No input fasta file for prefix $inputPrefix"
	fi

	if has_reference_index "$inputPrefix"; then 
		# First check if an index already exists
		warn "BWA index for input prefix ${inputPrefix} already exists. Could not build index"
		exitStatus=2
	else 	# Otherwise assume that input is a fasta-file and create an index
		info "Building BWA index from file $inputFasta"
		
		if [ -z "$outputPrefix" ]; then outputPrefix="$(basename $inputPrefix)"; fi

		bwaIndexArgs=(index -p "$outputPrefix" "$inputFastaFile")

		exec_command $BWA_CMD "${bwaIndexArgs[@]}" "1>&2"

		exitStatus=$?
	fi

	return $exitStatus
}

find_reference_index () {
	find_files "$1" "$BWA_INDEX_REGEX" 'BWA Index-'

	return $?	
}

has_reference_index () {
	local indexPrefix="$1"

	local indexFiles=( $( find_reference_index "$indexPrefix" ) )
	local numIndexFiles="${#indexFiles[@]}"

	debug "Wants $BWA_INDEX_NUMFILES files. Found $numIndexFiles"

	if [ "$numIndexFiles" -eq "$BWA_INDEX_NUMFILES" ]; then 
		printf '%s\n' "${indexFiles[@]}"
		return 0
	else
		return 1
	fi
}

map_reads () {
	local indexPrefix="$1"
	local inputPrefix="$2"
	local outputPrefix="$3"
	local readGroupString="$4"

	if [[ $# -lt 2 || $# -gt 4 ]]; then
		croak "Argument error"
	fi

	# Collect input files based on prefix
	local inputFiles=( $(find_fastq "$inputPrefix" ) )
	local numInputFiles=${#inputFiles[@]}
	
	if [[ $numInputFiles -ge 1 && $numInputFiles -le 2 ]]; then 
		debug "inputFiles: [ ${inputFiles[*]} ]"
	else
		croak "Wrong number of input files ($numInputFiles) found for prefix $inputPrefix";
	fi

	# Collect index files	
	local indexFiles=( $( has_reference_index "$indexPrefix" ) )
	
	if [ "$?" -ne 0 ]; then
		warn "No valid index found for reference prefix $indexPrefix";
		return 1
	fi

	# Build exec string
	bwaExecArgs=("$BWA_MAPPING_ENGINE" -t "$NUM_THREADS")

	# Query consists of interleaved reads
	if $BWA_QUERY_INTERLEAVED; then
		bwaExecArgs=("${bwaExecArgs[@]}" -p)
	fi

	# Initialize temporary folder
	init_work_folder "${inputFiles[@]}" "${indexFiles[@]}"

	# Strip folder names from input files
	indexPrefix="$(basename "$indexPrefix")"
	inputPrefix="$(basename "$inputPrefix")"

	inputFiles=( $( strip_paths "${inputFiles[@]}" ) )

	# Use input prefix if no output prefix was provided	
	if [ -z "$outputPrefix" ]; then
		outputPrefix="$inputPrefix"
	fi

	local outputFolder="$(dirname "$outputPrefix")"
	outputPrefix="$(basename "$outputPrefix")"

	# If no read group string was added, then add the minimum required information,
	# based on the input prefix
	if [ -z "$readGroupString" ]; then
		readGroupString="@RG\tID:${inputPrefix}\tSM:${outputPrefix}"
	fi

	# Do this or tabs are escaped incorrectly which messes up bwa execution
	#readGroupString=$(echo -e "$readGroupString" | sed 's/\t/\\t/g' )

	case $BWA_MAPPING_ENGINE in

	mem)
	    bwaExecArgs=("${bwaExecArgs[@]}" \
	                 -R "'$readGroupString'" \
	                 "$indexPrefix" \
	                 "${inputFiles[@]}" \
	    )
	    ;;
	aln)
	    croak "bwa backtrack mapping (aln) not implemented"
	    ;;
	bwasw)
	    croak "bwasw mapping not implemented"
	    ;;
	*)
	    croak "Unkown bwa mapping engine $BWA_MAPPING_ENGINE"
	esac

	# Execute read mapping
	exec_command $BWA_CMD "${bwaExecArgs[@]}" ">out/${outputPrefix}.sam"

	local exitStatus=$?

	# Wrap up and transfer output to the specified output folder
	close_work_folder "$outputFolder"

	return $exitStatus
}

# Report loading of module
debug "Plugin $(basename "${BASH_SOURCE[0]}") loaded"
debug "BWA version: $BWA_VERSION"
debug "BWA Index Regex: $BWA_INDEX_REGEX"
debug "BWA Number of expected index files: $BWA_INDEX_NUMFILES"
