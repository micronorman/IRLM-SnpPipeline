#!/usr/bin/env bash

#
#    This script is intended as a preliminary step for processing 
#    sequencing reads before they are ready to be mapped.
#    This usually involves removing traces of adapter fragments
#    as well as removing 'low-quality segments' from 5'- and 3'-ends.
#
#    Copyright (C) 2014 - 2016, Anders Norman
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

# Import common library functions
if [ ! -d "$SANDBOX_HOME/common" ]; then
	echo "Common library files not found" 1>&2
	exit 1
else
	for script in "$SANDBOX_HOME"/common/*.sh; do
		. $script
	done
fi

# Load plugins
. $SANDBOX_HOME/plugins/reads.sh
. $SANDBOX_HOME/plugins/bwa.sh
. $SANDBOX_HOME/plugins/samtools.sh
. $SANDBOX_HOME/plugins/picard.sh

# Hardcoded defaults for options with numeric values
defaultOutputFolder='.'
defaultOutputFormat='bam'
defaultNumThreads=$(( $NUM_THREADS / 2))

usage() {
	echo -e "Usage: $(basename $0) [ options ] reference_prefix query_prefix [ output_prefix ]\n\n"      \
		"       --help                 This help message\n"                                            \
		"       --sampleName <string>  This will add the given sample name to the Read Group (SM field).\n" \
		"       --outFormat [bam|cram] Choose bam or cram as the final output format [${defaultOutputFormat}].\n" \
		"       --minMapQ <integer>    The minimum acceptable mapping quality for keeping alignments [${READS_MINMAPQ}].\n" \
		"       --minBAQ <integer>     The minimum acceptable base quality [${READS_MINBAQ}].\n" \
		"       --noFixMates           Do not attempt to improve alignment by fixing mate-pair information.\n" \
		"       --noStats              Skip mapping statitics step. Speeds up mapping.\n" \
		"       --noMarkDupes          Do not mark likely PCR duplciate reads. Speeds up mapping time\n" \
		"       --quick                Only run the mapping step (equivalent to --noFixMates --noStats --noMarkDupes --noKeepUnmapped).\n"\
		"       --cleanSAM             Only keep mapped non-duplicate reads in the final BAM file\n" \
		"       --noKeepUnmapped       Do not store unmapped reads in seperate FASTQ-files if cleanSAM is set\n" \
		"       --plotStats            Make plots of mapping statistics. Will not run if --noStats is set.\n" \
		"       --verbose              Verbose output. Read all the things!\n" \
		"       --quiet                Suppress all output except errors.\n" \
		"       --log                  Redirect output to logfile (e.g. logfile.txt).\n" \
		"       --numThreads <n>       Number of threads to use for mapping [$defaultNumThreads].\n" \
		"       --dryRun               Don't actually run anything, but print all commands to be executed.\n" \
		"       --keepTemps            Do not delete the temporary folders that were created during the run.\n" \
		
	1>&2
	exit 1;
}

parse_prefixes () {
	# Locate existing files from the reference prefix and list all found input files
	local referencePrefix="$1"
	local queryPrefix="$2"

	local referenceFastaFile=
	local referenceFastaIndex=
	local referenceIndexFiles=()
	local referenceSequenceDictionary=
	local queryFiles=()
	local inputFileList=()

	if is_fasta $referencePrefix; then
		referenceFastaFile="$referencePrefix"
	else
		referenceFastaFile="$( find_fasta "$referencePrefix" )" 
	fi
	
	if [ -z "$referenceFastaFile" ]; then
		warn "No reference sequence file was found for reference prefix ${referencePrefix}"
	else 
		info "Using file $referenceFastaFile as the reference sequence"

		referenceFastaIndex="$( find_files "$referencePrefix" ".fai" )"
		referenceSequenceDictionary="$( find_files "$referencePrefix" ".dict" )"
	fi
	
	referenceIndexFiles=( $( has_reference_index "$referencePrefix" ) )
	queryFiles=( $( find_fastq "$queryPrefix" ) )

	if [ "${#queryFiles[@]}" -eq 1 ]; then
		info "Mapping single-end reads"
		noFixMates=true
		noMarkDupes=true
	else
		info "Mapping paired-end reads"
	fi	

	inputFileList=( "$referenceFastaFile" "$referenceFastaIndex" "$referenceSequenceDictionary"\
	                "${referenceIndexFiles[@]}" "${queryFiles[@]}" \
	)

	for filePath in "${inputFileList[@]}"; do
		[ -n "$filePath" ] && echo "$filePath"
	done	
}

run_mapping_pipe () {
	debug "$*"

	local referencePrefix="$1"
	local queryPrefix="$2"
	local outputPrefix="$3"
	local outputFormat="$4"
	local RG_string="$5"

	local inputFileList=( $( parse_prefixes "$referencePrefix" "$queryPrefix" ) )

	init_work_folder "${inputFileList[@]}"
 
	queryPrefix="$( basename "$queryPrefix" )"
	referencePrefix="$( basename "$referencePrefix" )"

	outputFolder="$( dirname "$outputPrefix")"
	outputPrefix="$( basename "$outputPrefix" )"

	if ! has_reference_index "$referencePrefix" > /dev/null; then
		build_reference_index "$referencePrefix" > /dev/null 
	fi

	local nextOutFile="$outputPrefix.sam"
	local lastOutFile=
	
	map_reads "$referencePrefix" "$queryPrefix" "$outputPrefix" "$RG_string"

	if [ $? -eq 0 ]; then
		lastOutFile="$nextOutFile"
	else
		warn "Mapping produced errors. Cannot continue with query $queryPrefix"
		return 1;
	fi

	# Fix mate-pair information
	if ! $noFixMates; then
		nextOutFile="$( relabel_filename "$lastOutFile" "fixedMates" )"

		fix_mates "$lastOutFile" "$nextOutFile"

		if [ $? -eq 0 ]; then
			rm "$lastOutFile"

			lastOutFile="$nextOutFile"
		else
			warn "Something went wrong when trying to fix mate-pair information"
		fi
	fi

	# Perform quick coordinate sorting
	local oldSamFmt=$SAMTOOLS_OUTPUT_FORMAT
	local oldSamLvl=$SAMTOOLS_COMPRESSION_LEVEL

	SAMTOOLS_OUTPUT_FORMAT='bam'
	SAMTOOLS_COMPRESSION_LEVEL=0

	nextOutFile="$( relabel_filename "$lastOutFile" "coordinateSorted" )"
	nextOutFile="$( basename "$nextOutFile" .sam )"
	
	coordinatesort_sam "$lastOutFile" "$nextOutFile"

	if [ $? -eq 0 ]; then
		rm "$lastOutFile"

		lastOutFile="$nextOutFile.bam"
		SAMTOOLS_OUTPUT_FORMAT=$oldSamFmt
		SAMTOOLS_COMPRESSION_LEVEL=$oldSamLvl
	else
		croak "Failed to perform coordinate sorting of alignments of $queryPrefix to $referencePrefix"
	fi

	if ! $noMarkDupes; then
		nextOutFile="$( relabel_filename "$lastOutFile" "markedDupes" )"
	
		mark_duplicates "$lastOutFile" "$nextOutFile"

		if [ $? -eq 0 ]; then
			rm "$lastOutFile"

			lastOutFile="$nextOutFile"
		else
			warn "Something went wrong when trying to mark duplicates"
		fi
	fi

	#if ! $noRealignIndels; then
	#	nextOutFile="$( relabel_filename "$lastOutFile" "gatkRealigned" )"
	#
	#	realign_indels "$lastOutFile" "$referencePrefix" "$nextOutFile"
	#
	#	rm "$lastOutFile"
	#
	#	lastOutFile="$nextOutFile"
	#fi

	if ! $noStats; then
		local statsOutFile="${outputPrefix}.mapping-stats.txt"

		generate_mapping_stats "$lastOutFile" "$referencePrefix" "out/$statsOutFile"

		if $plotStats; then
			stats_plot_bamstats "out/$statsOutFile" "out/stats/${outputPrefix}.mapping"
		fi
 	fi

	nextOutFile="out/$outputPrefix"
		
	if ! $noKeepUnmapped && [ "$SAMTOOLS_EXCLUDE_BITS" == "0xF04" ]; then
		store_unmapped_reads "$lastOutFile" "out/$outputPrefix"
	fi

	# Final step. Filters unmapped and low-quality reads
	finalize_sam "$lastOutFile" "$nextOutFile" "$referencePrefix" "$outputFormat"

	#mv "$lastOutFile" "out/$outputPrefix.bam"

	close_work_folder "$outputFolder"
}

# Check user input
if [[ $# -lt 1 || "$1" == "--help" ]]; then usage; fi

noFixMates=false
noMarkDupes=false
noStats=false
noKeepUnmapped=false
plotStats=false
oldNumThreads=$NUM_THREADS

NUM_THREADS="$defaultNumThreads"

while [ $# -gt 0 ]; do

	case "$1" in
		--sampleName)
			check_pathsafe "$2" "$1"
			sampleName="$2"; shift 2;;
		--outFormat)
			outFormat="$2"; shift 2;;
		--minMapQ)
			check_posinteger "$2" "$1"
			READS_MINMAPQ="$2"; shift 2;;
#		--noRealignIndels)
#			noRealignIndels=true; shift;;
		--noFixMates)
			noFixMates=true; shift;;
		--noMarkDupes)
			noMarkDupes=true; shift;;
		--noKeepUnmapped)
		    noKeepUnmapped=true; shift;;
		--noStats)
			noStats=true; shift;;
		--plotStats)
			plotStats=true; shift;;
		--cleanSAM)
			SAMTOOLS_EXCLUDE_BITS="0xF04"; shift;;
		--quick)
			noFixMates=true; 
			noMarkDupes=true;
			noStats=true;
			noKeepUnmapped=true;
			shift;;
		--numThreads)
			check_posinteger "$2" "$1"
			NUM_THREADS="$2"; shift 2;;
		--dryRun)
			NO_EXEC=true; shift;;
		--keepTemps)
			KEEP_TEMPS=true; shift;;
		--verbose)
			set_log_level 3; shift;;
		--log)
			check_pathsafe "$2" "$1"
			log_file "$2"; shift 2;;
		--quiet)
			set_log_level 0; MUTE_EXEC=true; MUTE_CMD=true; shift;;
		--debug)
			set_log_level 4; KEEP_TEMPS=true; shift;;
		--)
			shift; break;;
		--*)
			echo "Unrecognized option $1" 1>&2; usage;;
		-*)
			echo "Unrecognized option $1" 1>&2; usage;;
	*)
		break;
	esac
done

# Check user input and set input variables
if [[ $# -lt 2 || $# -gt 3 ]]; then echo "Wrong number of input arguments" 1>&2; usage; fi

referencePrefix="$1"
inputPrefix="$2"
outputPrefix="$3"
platform="ILLUMINA"

if [ -z "$outputPrefix" ]; then
	outputPrefix="$defaultOutputFolder/$(basename "$inputPrefix")"
fi

if [ -z "$outFormat" ]; then
	outFormat="${defaultOutputFormat}"
fi

if [ -n "$sampleName" ]; then
	readGroupHeaderString="@RG\\tID:$(basename ${inputPrefix})\\tSM:${sampleName}\\tPL:${platform}"
fi


run_mapping_pipe "$referencePrefix" "$inputPrefix" "$outputPrefix" "$outFormat" "$readGroupHeaderString"

NUM_THREADS=$oldNumThreads

exit 0
