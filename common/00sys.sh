#!/usr/bin/env sh

# The operating system runnning the script.
OPSYS="$(uname -s)"

if [ "$OPSYS" = 'Linux' ]; then
	REALPATH="$(type -p readlink)"
	REALPATH="$REALPATH -f"
elif [ "$OPSYS" = 'Darwin' ]; then
	REALPATH="$(type -p grealpath)"
fi

# The maximum number of threads allowed for any executed subshell command
if [ -z "$NUM_THREADS" ]; then 
	if [ "$OPSYS" = "Darwin" ]; then
		NUM_THREADS="$( sysctl hw.ncpu | cut -f 2 -d' ' )"
	elif [ "$OPSYS" = "Linux" ]; then
		NUM_THREADS="$( grep -c processor /proc/cpuinfo )"
	fi
fi

# Available system memory
if [ "$OPSYS" = 'Linux' ]; then
	SYSMEM="$(cat /proc/meminfo | grep MemTotal | sed 's/MemTotal:[ ^I]*\([0-9]*\) kB$/\1/')"
elif [ "$OPSYS" = 'Darwin' ]; then
	SYSMEM="8589935"
fi

# Native script paths
SCRIPT_PATH="$( $REALPATH "$0" )"
SCRIPT_NAME="$(basename "$0")"
SCRIPT_HOME="$(dirname "$0")"
RUN_DIR="$(pwd)"

# Set the debug flag (disabled by default)
if [ -z "$DEBUG" ]; then
	DEBUG=false
fi

if [ -z "$VERBOSE" ]; then
	VERBOSE=false
fi

if $DEBUG; then
	VERBOSE=true
fi

HAS_LOG=false

## Regular expressions to verify variable types

# Numeric (pos/neg) float with optional exponent
REGEX_NUMERIC="^[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?$"

# Positive integer
REGEX_POSINTEGER="^[0-9]+$"

# Pathsafe
REGEX_PATHSAFE="^[-.0-9A-Z_a-z]+$"

# Nucleotides
REGEX_NT="^[ATCGNatcgn]+$"


croak () {
	# Gracefully exit with an error message
	local msg=${1:-"Unkown error"} >&2
	local status=${2:-1}

	if $DEBUG; then
		msg="[ ${SCRIPT_NAME}::${FUNCNAME[1]} ] $msg"
	fi

	# Use the log function when possibly. Default to echo
	$HAS_LOG && log "$msg" 0 "ERROR" || echo "ERROR: $msg" >&2

	exit $status	
}

warn () {
	# A non-fatal warning sent to stderr
	[ -n $HAS_LOG ] && log "$1" 1 "WARNING"|| echo "WARNING: $1" >&2
}

message () {
	# A bare message for the user sent to stderr. Will be muted in quiet mode
	[ -n $HAS_LOG ] && log "$1" 2 || echo "$1" >&2
}

info () {
	# Verbose information sent to stderr. Will only be displayed in verbose mode
	if ! $VERBOSE; then
		return 0
	elif $HAS_LOG; then
		log "$1" 3
	else
		echo "$1" >&2
	fi
}

debug () {
	if ! $DEBUG;then
		return 0
	fi

	local msg="[ ${SCRIPT_NAME}::${FUNCNAME[1]} ] $1"

	if $HAS_LOG; then
		log "$msg" 4
	else
		echo "$msg" >&2
	fi
}

join () {
	# Join an array into a text string with a specified delimiter
	local IFS="$1"; shift; echo "$*";
}

exec_path () {
	$REALPATH $(type -p "$1") 2> /dev/null
}

strip_paths () {
	for pathName in "$@"; do
		basename "$pathName"
	done
}

get_extension () {
	local fileName="$(basename "$1")"

	if [ -z "$fileName" ]; then return 1; fi

	echo "${fileName##*.}"
}

relabel_filename () {
	local fileName="$1"
	local label="$2"

	
	local fileExt="$( get_extension "$fileName" )"

	if [ -z "$fileExt" ]; then
		echo "${fileName}.${label}"
	else
		echo "${fileName}" | sed "s/${fileExt}/${label}.${fileExt}/"
	fi
}

real_paths () {
	for fileName in "$@"; do
		$REALPATH "$fileName"
	done
}

find_files () {
	debug "$*"
	local prefix="$1"
	local suffix="$2"
	local label="$3"

	local inFolder="$(dirname "$prefix")"
	prefix="$(basename "$prefix")"

	fileList=( $( find "$inFolder/" -maxdepth 1  -print 2> /dev/null | sort | grep -E "\/$prefix.*$suffix$" 2>/dev/null ) )

	if [ "${#fileList[@]}" -gt 0 ]; then
		echo "${fileList[*]}"
	else
		debug "No files found in folder $inFolder for prefix $prefix and suffix $suffix"
		return 1
	fi

	return 0
}

check_numeric () {
	[[ "$1" =~ $REGEX_NUMERIC ]] || croak "Value for $2 ($1) is invalid, must be numerical"
}

check_pathsafe () {
	[[ "$1" =~ $REGEX_PATHSAFE ]] || croak "Value for $2 ($1) contains one or more non-pathsafe characters"
}

check_posinteger () {
	[[ "$1" =~ $REGEX_POSINTEGER ]] || croak "Value for $2 ($1) is invalid, must be a positive integer"
}

check_exists () {
	# Check if a file exists. Croaks if not. Overridden in NO_EXEC mode
	# Second argument should indicate the type of file being looked for (e.g. FASTA)
	if [ -z $1 ]; then
		croak "Empty filename given with which to check existence of $2"
	elif [ ! -e $1 ];then
		if $NO_EXEC; then return 0; fi
		croak "${2:-"Unkown"} file $1 not found"
	fi

}

debug "Debug mode: on"
debug "Verbose mode: on"
debug "Operating system: $OPSYS" 
debug "Available memory (kB): $SYSMEM"
debug "Number of CPUs: $NUM_THREADS"
debug "Script real path: $SCRIPT_PATH"
debug "Script name: $SCRIPT_NAME"
debug "Script folder: $SCRIPT_HOME"
debug "Run dir: $RUN_DIR"
debug "Loaded system module $(basename ${BASH_SOURCE[0]})"
