#!/usr/bin/env sh

# Keep track of all temporary files created while the script was running
NUM_TEMPS=0
TEMP_DIRS_LIST=()
TEMP_INDEX=
CURR_WORK_DIR="$RUN_DIR"
PREV_WORK_DIR=
NO_CLOBBER=true

# Do not delete temporary directories at exit
KEEP_TEMPS=false

# Command for creating temporary directory.
MKTEMP="$(type -p mktemp) -dt"

# Resolve different behaviour from mktemp under Darwin and Linux
if [ "$OPSYS" = 'Darwin' ]; then
	MKTEMP_TEMPLATE="tmp.$USER.$(echo "$SCRIPT_NAME" | tr . _)"
elif [ "$OPSYS" = 'Linux' ]; then
	MKTEMP_TEMPLATE="tmp.$USER.$(echo "$SCRIPT_NAME" | tr . _ )".XXXXXX
else
	croak "Unsupported operating system $OPSYS"
fi

# Setup trap to remove all temporary folders at exit
trap "debug 'Program exit' && cleanup" EXIT

make_temp () {
	# A unique temporary directory will be created based on the script and user names
	# All temporary directories will be removed at program exit unless the KEEP_TEMPS 
	# flag has been set
	local tempDir="$($MKTEMP $MKTEMP_TEMPLATE)"
	local exitStatus=$?

	if [ $exitStatus -ne 0 ]; then
		warn "Failed to create temporary directory. Mktemp exit status: $exitStatus"
		return $exitStatus
	else
		# Add new temporary directory to internal list. These will all be removed during cleanup
		TEMP_DIRS_LIST+=("$tempDir")
		let NUM_TEMPS++

		debug "Created temporary directory $tempDir"
		debug "Number of temporary directories: $NUM_TEMPS"
	fi
	
	return 0
}

remove_temp () {
	# First, check that there are any temp directories left
	if [ "$NUM_TEMPS" -gt 0 ]; then
		let NUM_TEMPS--
		tempDir="${TEMP_DIRS_LIST[$NUM_TEMPS]}"

		if ! $KEEP_TEMPS; then

			# As a precaution, only directories beginning with 'tmp.' can be
			# deleted by this routine
			tmpDirName="$(basename "$tempDir")"

			if [ ! ${tmpDirName:0:4} == "tmp." ]; then
				warn "$SCRIPT_NAME tried to remove a non-temp dir"
			else
				# Here rm is called on the temporary directory
				rm -rf "$tempDir"
			
				if [ $? -eq 0 ];then
					debug "Removed temporary directory $tempDir"
				else
					 warn "Failed to delete temporary directory $tempDir"
				fi
			fi
		fi

		unset TEMP_DIRS_LIST[$NUM_TEMPS]
	else 
		warn "remove_temp was called with no temporary directories left"
	fi

	debug "Number of temporary directories left: $NUM_TEMPS"
}

enter_temp () {

	# Set or increase the index of the temporary directory list
	if [ -z "$TEMP_INDEX" ]; then
		TEMP_INDEX=0
	else 
		let TEMP_INDEX++
	fi

	debug "TEMP_INDEX: $TEMP_INDEX"
	
	# Change to new temporary directory or complain about not being able to
	if [ "$TEMP_INDEX" -lt "$NUM_TEMPS" ]; then

		local nextTempDir="${TEMP_DIRS_LIST[$TEMP_INDEX]}"

		if [ ! -d "$nextTempDir" ]; then
			croak "Temporary directory $nextTempDir no longer exist"
		fi
	
		PREV_WORK_DIR=$(pwd)
		CURR_WORK_DIR=$nextTempDir

		cd $nextTempDir
		
		debug "Entering $nextTempDir. PREV_WORK_DIR is now $PREV_WORK_DIR"
	else
		warn "No more temporary directories to enter. Staying put"
		let TEMP_INDEX--
	fi

	debug "CURR_WORK_DIR: $CURR_WORK_DIR"
	debug "PREV_WORK_DIR: $PREV_WORK_DIR"
}

leave_temp () {

	if [ -n "$TEMP_INDEX" ]; then
	
		# Unset or decrease the temporary index counter
		if [ $TEMP_INDEX -le 0 ]; then
			unset TEMP_INDEX
		else 
			let TEMP_INDEX--
		fi		

		debug "TEMP_INDEX: $TEMP_INDEX"

		if [ -n "$TEMP_INDEX" ]; then
			prevTempDir="${TEMP_DIRS_LIST[$TEMP_INDEX]}"
			

			debug "leaving $CURR_WORK_DIR and entering $prevTempDir"

			CURR_WORK_DIR="$prevTempDir"

			if [ $TEMP_INDEX -eq 0 ]; then
				PREV_WORK_DIR="$RUN_DIR"
			else 
				PREV_WORK_DIR="${TEMP_DIRS_LIST[$((TEMP_INDEX-1))]}"
			fi
		else
			debug "This was the last of the temporary directories"

			CURR_WORK_DIR="$RUN_DIR"
			unset PREV_WORK_DIR
		fi
	else
		warn "Tried to leave a temp folder when not in a temp folder"
		CURR_WORK_DIR="$RUN_DIR"
		unset PREV_WORK_DIR
	fi

	
	# Finally, we can change directory
	if [ ! -d "$CURR_WORK_DIR" ]; then
		croak "Temporary directory $CURR_WORK_DIR no longer exists"
	else 
		cd "$CURR_WORK_DIR"
		debug "CURR_WORK_DIR: $CURR_WORK_DIR"
		debug "PREV_WORK_DIR: $PREV_WORK_DIR"
	fi
}

init_work_folder () {
	if $NO_EXEC; then return; fi

	debug "New work folder requested( $# files )"
	
	make_temp
	enter_temp

	# Link input files to current working directory
	for inputFile in "$@"; do
		realPath="$($REALPATH "$PREV_WORK_DIR/$inputFile")"
		symLink="$(basename "$inputFile")"

		# If full path is empty, assume that it was not a relative path
		# and attempt to link to the absolute path instead
		if [ -z "$realPath" ]; then
			realPath="$($REALPATH "$inputFile")"
		fi

		if [ ! -e "$symLink" ]; then
			ln -s "$realPath" "$symLink"
			
			if [ $? -eq 0 ]; then
				debug "Linked input file $inputFile to folder $CURR_WORK_DIR"
			else
				warn "Failed to link $inputFile to folder $CURR_WORK_DIR"
				return 1
			fi 
		else
			warn "Attempted to link '$inputFile' to existing file or folder '$symLink'. Skipping"  
		fi
	done

	# Create a folder to collect output
	mkdir -p 'out'
}

add_to_work_folder () {
	if $NO_EXEC; then return; fi

	# Check if we are in a temporary working folder 
	if [ -z "$PREV_WORK_DIR" ]; then
		warn "We don't apear to be in a working folder"
		return 1
	fi

	debug "Adding $# files or directories to $CURR_WORK_DIR"

	# Link additional input files to current working directory
	for inputFile in "$@"; do
		realPath="$($REALPATH "$PREV_WORK_DIR/$inputFile")"
		symLink="$(basename "$inputFile")"

		if [ ! -e "$symLink" ]; then
			ln -s "$realPath" "$symLink" 2>/dev/null
			
			if [ $? -eq 0 ]; then
				debug "Linked input file $inputFile to folder $CURR_WORK_DIR"
			else
				debug "Failed to link $inputFile to folder $CURR_WORK_DIR (ln exit status: $?)"
				return 1
			fi   
		fi
	done
}

close_work_folder () {
	if $NO_EXEC; then return; fi

	debug "Closing work folder $CURR_WORK_DIR"

	local outputFolder="$1"
	shift

	if [ -z "$outputFolder" ]; then
		outputFolder='.'
	fi

	mkdir -p "$PREV_WORK_DIR/$outputFolder"

	# Check that output folder has contents	
	if [ -d 'out' ]; then
		cd 'out'

		debug "Output files: $( ls -m )"

		if [ "$PREV_WORK_DIR" != "$RUN_DIR" ]; then
			
			# Since we are inside temporary folders use hard links instead of cp
			# which is quicker
			ln * "$PREV_WORK_DIR/$outputFolder" 2> /dev/null

			if [ $? -eq 0 ]; then
				debug "Hardlinked output files in $CURR_WORK_DIR/out to $PREV_WORK_DIR/$outputFolder"
			else	
				warn "Failed to link output files in $CURR_WORK_DIR/out to $PREV_WORK_DIR/$outputFolder"
			fi
		else	
			# copy (without clobber) output files to destination
			cp -nrp "." "$RUN_DIR/$outputFolder"

			if [ $? -eq 0 ]; then
				debug "Transferred output files in $CURR_WORK_DIR/out to $RUN_DIR/$outputFolder"
			else
				warn "Failed to transfer output files in $CURR_WORK_DIR/out to $RUN_DIR/$outputFolder"
			fi
		fi

		cd ..

	else
		warn "Working folder $CURR_WORK_DIR contained no output"		
	fi

	leave_temp
	remove_temp	
}

cleanup () {
	# Removes all temporary directories unless KEEP_TEMPS is set
	debug "Cleanup begin"

	if [ "${#TEMP_DIRS_LIST[@]}" -gt 0 ]; then 
		debug "Temporary directories left: ${TEMP_DIRS_LIST[*]}"
	fi

	# Try not to interrupt the cleanup process even if Ctrl-C is hit
	trap 'warn "Cleanup in progress. No stopping now!"' INT TERM

	if ! $KEEP_TEMPS; then
		while [ "${#TEMP_DIRS_LIST[@]}" -gt 0 ]; do
			remove_temp
		done
	else
		debug "Keeping temporary folders. Nothing to be done"
	fi

	if [ -d "$RUN_DIR" ]; then
		cd "$RUN_DIR"
	else 
		warn "Directory $RUN_DIR no longer exists and could not be deleted"
		cd "$HOME"
	fi

	debug "Cleanup after $SCRIPT_NAME completed @ $(date)"

	trap INT TERM
}

debug "System library $(basename "${BASH_SOURCE[0]}") loaded"
debug "Base mktemp template: $MKTEMP_TEMPLATE" 
