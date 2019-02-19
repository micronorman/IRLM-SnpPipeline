#!/usr/bin/env bash

JAVA_HEAP_MEM='4G'

# Execution of subshells

if [ -z "$NO_EXEC" ]; then NO_EXEC=false; fi
if [ -z "$MUTE_CMD" ]; then MUTE_CMD=false; fi

# Except in debug mode, output of executed commands is suppressed
if [ -z "$MUTE_EXEC" ]; then
	if $DEBUG; then
		MUTE_EXEC=false
		MUTE_CMD=false
	else 
		MUTE_EXEC=true
	fi
fi

# Global nice level. Executed command will be run with lower nice-value by default
if [ -z "$NICENESS" ]; then
	NICENESS=10
fi

if [ -z "$JAVA_BIN" ]; then
	JAVA_BIN="$( exec_path java )"
fi

echo_command () {
	if ! $MUTE_CMD; then
		$HAS_LOG && log "$*" 2 || echo "$*" >&2
	fi
}

exec_command () {
	local cmdName="$(basename "$1")"
	local exitStatus=0

	echo_command "$@"

	if ! $NO_EXEC; then

		# Where to direct command output
		if $MUTE_EXEC; then
			debug "Redirecting $cmdName output to /dev/null"
			exec 3>/dev/null
		else
			debug "Redirecting $cmdName output to stderr"
			exec 3>&2
		fi

		# Trap Ctrl-C to ensure proper cleanup
		trap "exec 3>&- && croak 'Program $cmdName stopped by user' 130" INT TERM

		# Execute the command. NOTE: Apparently, executing with eval is unsafe...
		eval "nice" "-n" "$NICENESS" "$@" 2>&3

		exitStatus=$?
		
		echo >&3
		exec 3>&-
		
		debug "$cmdName exit status: $exitStatus"
	fi

	
	if [ $exitStatus -ne 0 ]; then
		warn "$cmdName finished with errors (exit status: $exitStatus)"
	fi


	return $exitStatus
}

exec_java_command () {
	local jarFile="$1"
	
	check_exists "$jarFile" "Java archive"
	shift;

	local javaCmdArgs=( -Xmx${JAVA_HEAP_MEM} -jar $jarFile )

	exec_command $JAVA_BIN "${javaCmdArgs[@]}" "$@"
	exitStatus=$?

	return $exitStatus 
}

debug "System library $(basename "${BASH_SOURCE[0]}") loaded"
debug "Niceness level: $NICENESS"
debug "Java binary: $JAVA_BIN"

if $MUTE_CMD; then debug "Mute commands"; else debug "Output commands"; fi
if $MUTE_EXEC; then debug "Mute command output"; else debug "Output command output"; fi


