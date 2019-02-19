#!/usr/bin/env bash

#    pre_trim_reads - A general purpose script for the pre-proccessing of 
#    Illumina-based sequenceing reads. It's really just a wrapper for 
#    Trimmomatic (see <http://www.usadellab.org/cms/?page=trimmomatic>)
#
#    This script is intended as a preliminary step for processing 
#    sequencing reads before they are ready to be mapped.
#    This usually involves removing traces of adapter fragments
#     as well as removing 'low-quality segments' from 
#
#    Copyright (C) 2014-2018, Anders Norman
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
. $SANDBOX_HOME/plugins/trimmomatic.sh

defaultInputFolder='.'
defaultOutputFolder='.'

usage() {
	echo -e "Usage: $(basename $0) [ options ] input_prefix sample_prefix\n\n" \
		"       --help                      This help message\n" \
		"       --noClipAdapters            Don't remove adapter fragments (not recommended) [disabled]\n" \
		"       --adaptersPrefix            Choose beteween TruSeq2,TruSeq3,NexteraPE & All [ All ]\n" \
		"       --clipAdaptSeedMismatches   The maximum mismatch count to the adapter sequence [$TRIMMOMATIC_CLIPADAPT_SEEDMM]\n" \
		"       --clipAdaptSimpleThreshold  Accuracy of 'simple' matches to the adapter [$TRIMMOMATIC_CLIPADAPT_SIMPLETHR]\n"    \
		"       --clipAdaptPalinThreshold   Accuracy of match between 'adapter ligated reads' [$TRIMMOMATIC_CLIPADAPT_PALINTHR]\n"   \
		"       --endTrimLeadQual           Trim leading bases (5') below the indicated phred-score [$TRIMMOMATIC_ENDTRIM_LEADQUAL]\n" \
		"       --endTrimTrailQual          Trim trailing bases (3') below the indicated phred-score [$TRIMMOMATIC_ENDTRIM_TRAILQUAL]\n" \
		"       --filterMinLen              Filter reads shorter than this length (after trimming) [$TRIMMOMATIC_FILTER_MINLEN]\n" \
		"       --noConvertPhred64          Prevents conversion of legacy base64 phred-scores into base33 (not recommended) [ disabled ]\n" \
		"       --numThreads                Number of threads to be used for trimming [$NUM_THREADS]\n" \
		"       --keepUnpaired              Retain files of unpaired reads generated by Trimmomatic [ disabled ]\n" \
		"       --verbose                   Verbose output\n" \
		"       --log                       Redirect (verbose) output to a logfile (e.g. logfile.txt)\n" \
		"       --dryRun                    Don't actually run anything, but print all commands to be executed\n" \
		"       --keepTemps                 Do not delete the temporary folders that are created during the run\n" \
		
	1>&2
	exit 1;
}

# Check user input
if [[ $# -lt 1 || "$1" == "--help" ]]; then usage; fi

while [ $# -gt 0 ]; do
	case "$1" in
		--endTrimLeadQual)
			check_posinteger "$2"
			TRIMMOMATIC_ENDTRIM_LEADQUAL="$2"; shift 2;;
		--endTrimTrailQual)
			check_posinteger "$2"
			TRIMMOMATIC_ENDTRIM_TRAILQUAL="$2"; shift 2;;
		--numThreads)
			NUM_THREADS="$2"; shift 2;;
		--filterMinLen)
			check_posinteger "$2"
			TRIMMOMATIC_FILTER_MINLEN="$2"; shift 2;;
		--adaptersPrefix)
			adaptersPrefix="$2"; shift 2;;
		--noClipAdapters)
			TRIMMOMATIC_CLIPADAPT=false; shift;;
		--clipAdaptSeedMismatches)
			check_posinteger "$2"
			TRIMMOMATIC_CLIPADAPT_SEEDMM="$2"; shift 2;;
		--clipAdaptSimpleThreshold)
			check_posinteger "$2"
			TRIMMOMATIC_CLIPADAPT_SIMPLETHR="$2"; shift 2;;
		--clipAdaptPalinThreshold)
			check_posinteger "$2"
			TRIMMOMATIC_CLIPADAPT_PALINTHR="$2"; shift 2;;
		--noConvertPhred64)
			TRIMMOMATIC_TOPHRED33=false; shift;;
		--dryRun)
			NO_EXEC=true && set_log_level 2; shift;;
		--keepUnpaired)
			TRIMMOMATIC_KEEP_UNPAIRED=true; shift;;
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
			set_log_level 4; shift;;
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


#if [[ $# -lt 1 || $# -gt 3 ]]; then echo "Wrong number of input arguments" 1>&2; usage; fi

inputPrefix=$1
outputPrefix=$2

inputFiles=( $( find_fastq "$inputPrefix" ) )


if [ "${#inputFiles[@]}" -eq 1 ]; then
	info "Processing single-end reads"
	pre_process_se_reads "$inputPrefix" "$outputPrefix"
else
	info "Processing paired-end reads"
	pre_process_pe_reads "$inputPrefix" "$outputPrefix" "$adaptersPrefix"
fi	

