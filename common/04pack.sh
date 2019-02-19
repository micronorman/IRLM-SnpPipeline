#!/usr/bin/env sh

# Paths to pigz or gzip
PIGZ="$( exec_path pigz )"

# Use the following compression level (0-9) when compressing
GZIP_COMPRESS_LEVEL=6
GZIP_KEEP_SOURCE=false

# Use the multi-threaded program pigz as an alternative to gzip if installed
if [ -n "$PIGZ" ]; then
	GZIP_CMD="$PIGZ -p $NUM_THREADS"
else
	GZIP_CMD="$( exec_path gzip )"
fi

is_gzipped () {
	# Checks whether input argument is a gzipped file
	if [ -z "$1" ]; then croak "No input defined"; fi
	if [ ! -e "$1" ]; then croak "Input file $1 not found"; fi

	if $(file --mime-type "$($REALPATH "$1")" | grep -q gzip$); then
		return 0
	else
		return 1
	fi
}

gzip_unpack () {
	local inFile="$1"
	local outFile="$2"

	if ! is_gzipped "$inFile"; then croak "$inFile is not zipped"; fi

	if [ -z "$outFile" ]; then
		outFile="${inFile%.gz}"
	fi

	if [ "$inFile" == "$outFile" ]; then
		croak "Cannot unzip file $inFile unto itself"
	fi

	exec_command $GZIP_CMD -cdfq "$inFile" > "$outFile"
}

gzip_pack () {
	local inFile="$1"
	local outFile="$2"

	if is_gzipped "$inFile"; then croak "$inFile is already zipped"; fi

	if [ -z "$outFile" ]; then
		outFile="${inFile}.gz"
	fi

	if [ "$inFile" == "$outFile" ]; then
		croak "Cannot zip file $inFile unto itself"
	fi

	exec_command $GZIP_CMD -cfq "-${GZIP_COMPRESS_LEVEL}" "$inFile" > "$outFile"

	if ! $GZIP_KEEP_SOURCE; then
		rm $inFile;
	fi
}

debug "System library $(basename "${BASH_SOURCE[0]}") loaded"
debug "Compression command used: $GZIP_CMD"
