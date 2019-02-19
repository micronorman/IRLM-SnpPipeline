#!/usr/bin/env sh

# Logging levels
HAS_LOG=true

# The default logging level. Errors and warnings only
if $DEBUG; then 
	LOG_LEVEL=4
elif [ -z "$LOG_LEVEL" ]; then
	LOG_LEVEL=2 
fi

# By default, redirect the log stream to stderr
exec 4>&2

log () {
	local msg="$1"
	local level="$2"
	local label="$3"

	if [ -z "$msg" ]; then return; fi
	
	local verboseLevel=${level:-$LOG_LEVEL}
	
	if [ "$verboseLevel" -gt 4 ]; then
		croak "Illegal log level $verboseLevel"
	fi

	if [ "$LOG_LEVEL" -ge "$verboseLevel" ]; then
		if [ -n "$label" ]; then
			msg="$label: $msg"
		fi

		#local scriptTag="[ ${SCRIPT_NAME} ]:"

		#if $DEBUG; then
		#	scriptTag="[ ${scriptTag}::${FUNCNAME[2]} ]:"
		#fi

		#if [ -n "$verboseLevel" -ne 3 ]; then
		#	msg="$label [ $scriptTag ]: $msg"
		#fi
	
		echo -e "$msg" >&4
	fi

}

set_log_level () {
	local newLevel="$1"

	if [[ "$newLevel" -lt 0 || "$newLevel" -gt 4 ]]; then
		croak "User attempted to set illegal verbosity level ($newLevel)"
	fi

	case $newLevel in

	[01])
	# Log levels 0-1: Suppress command statement and output of executed commands
		DEBUG=false
		VERBOSE=false
		MUTE_EXEC=true
		MUTE_CMD=true
		;;
	2)
	# Log level 2: Show only Errors, warnings, and messages
		DEBUG=false
		VERBOSE=false
		MUTE_CMD=false
		MUTE_EXEC=true
		;;
	3)
	# Like level 3: But "verbose" information sent to output as well
		DEBUG=false
		VERBOSE=true
		MUTE_CMD=false
		MUTE_EXEC=false
		;;
	4)
	# Debug mode. All debug messages are sent to output. Executed command output is no longer hidden
		DEBUG=true
		VERBOSE=true
		MUTE_CMD=false
		MUTE_EXEC=false
		;;
	*)
		croak "Attempted to set illegal logging level ($newLevel) use levels 0-4"

	esac
	
	LOG_LEVEL="$newLevel"
	message "New log level: $newLevel"
}

log_header () {
	local msg="$1"

	border="$(printf '=%.0s' {1..80})"
	logString="$border\n$msg\n$border\n"

	log "$logString" 3
}

log_timed_event () {
	local msg="$1"
	local timeStamp="$(date)"

	log_header "$msg @ $timeStamp"
}

log_file () {
	local log_file="$( $REALPATH "$1" )"

	touch "$log_file" || croak "Log file $log_file is not writable"
	
	log "[ $SCRIPT_NAME ]: Redirecting output to $log_file" 3
	
	exec 4>>$log_file
}


if [ -n "$LOG_FILE" ]; then
	log_file "$LOG_FILE"
fi	

debug "System library $(basename "${BASH_SOURCE[0]}") loaded"
debug "Operating system: $OPSYS" 
debug "Maximum number of process threads: $NUM_THREADS"
debug "Script real path: $SCRIPT_PATH"
debug "Script name: $SCRIPT_NAME"
debug "Script folder: $SCRIPT_HOME"
debug "Run dir: $RUN_DIR"
debug "Logging level: $LOG_LEVEL"
