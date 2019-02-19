# NGS pipeline for calling variants (SNPs) of *Mycobacterium tuberculosis* at IRLM Copenhagen, Denmark

## Author
[Anders Norman](https://twitter.com/anders_norman)

## Equipment 
We normally use an off-site computer cluster managed by The Danish Health Authority (Sunhedsdatastyrelsen) that runs on Linux (Red Hat Enterprise Linux Server release 7.6). However, the pipeline has also been tested on a MacBook Pro (2.3 Ghz, Retina, Mid 2012 edition with 8GB RAM) running Mac OS X 10.12 (Sierra). 

## Requirements
* Perl >= 5.12
* Bioperl >= 1.6
* Vcftools >= 0.1.14
* R >= 3.0.0
* R package [ape](https://cran.r-project.org/web/packages/ape/) (Analysis of Phylogenetics and Evolution)

## Core tools (must be pre-installed)
* bwa (v0.7.15)
* samtools (v1.5)
* bcftools (v1.5)
* picard tools (v2.12.0)
* trimmomatic (v0.38)

These can easily be installed using [HomeBrew](http://brew.sh/) (Mac OS X) or [LinuxBrew](http://linuxbrew.sh/) (Linux), although be advised that we have only tested the pipeline up until and including the versions mentioned above.

## Pipeline scripts (located in the 'scripts' subfolder)
* pre_process_reads.sh
* map_reads.sh
* make_mapping_report.pl
* call_variants.pl
* vcf2fasta.pl
* aln2distmat.R

## Example

### Pointing to the home directory
In order to work, the environmental variable SANDBOX_HOME must first be set to the directory in which these files are located, like this:
```bash	
% export SANDBOX_HOME=`pwd`
```

### Example data
As an example of how the pipeline works, we are going to focus on four samples from the large Danish Cluster-2 outbreak (aka. 1112-15). The data is located in the NCBI [sequence read archive (SRA)](https://www.ncbi.nlm.nih.gov/sra)), and can be downloaded using the tool `fastq-dump` that comes with the SRA Toolkit:

Run accession | Sample name
--------------|-------------
ERR1949967 | DKC2-PP1
ERR1949970 | DKC2-PP15
ERR1949984 | DKC2-T101
ERR1950078 | DKC2-T91

### Pre-proccesing raw reads
Raw reads are first processed using the script pre_process_reads.sh, which directs trimmomatic to trim
away any remaining Illumina NexteraXT adapter fragments, and trailing low-quality (Q-score < 3) bases, like this:

```bash
% fastq-dump --gzip --split-files -O example ERR1949967

% scripts/pre_process_reads.sh example/ERR1949967 example/ERR1949967.trim
```

### Mapping reads to the *Mycobacterium tuberculosis* reference genome
We then execute the script `map_reads.sh`, which directs the program bwa to map reads to the *Mycobacterium tuberculosis H37Rv* reference genome with accession no. NC_000962 (version 3), located in the and subsequently uses picard-tools to mark PCR and optical duplicate reads. Finally, samtools is used to sort the resulting bam-file.

```bash
% scripts/map_reads.sh --sampleName DKC2-PP1 reference/NC_000962 ./ERR1949970.trim example/DKC2-PP1/ERR1949970
```

Note that we create a folder with the sample name, and use the prefix corresponding to the accession number for the mapping-file output. In this way, if a sample contains multiple runs, we can keep track of them within the same sample folder. We should now have two files: `ERR1949967.mapping-stats.txt` and `ERR1949967.bam`

### Assessing quality of the sample
Now we use the script `make_mapping_report.pl` to analyze the genome coverage in the sample folder, like this:
```
% scripts/make_mapping_report.pl -T reference/NC_000962.fna.gz -R reference/NC_000962.default_non-repetetive.bed example/DKC2-PP1 > example/DKC2-PP1/DKC2-PP1.coverage-report.txt
```
We use the bed-file `NC_000962.default_non-repetetive.bed` to disregard coverage in all regions of H37Rv that are later filtered out in downstream analysis. In this way, poorly covered regions do not affect the parts of the genome that are actaully analyzed. The output should look like this:
```
version: 0.4.11
sample: DKC2-PP1
number_of_libraries: 1
library_names: ERR1949967
library_number_of_reads: 2664599
library_mapping_frequencies: 0.979701636156135
library_error_rates: 9.511607e-03
bases_mapped: 286507282
positions_covered: 4005012
average_coverage: 71.5371844079368
stdev_coverage: 19.7067702269032
coverage_outside_2sd: 0.0552283938109468
minimum_coverage: 0
maximum_coverage: 285
coverage_rating: 3
coverage_thresholds:
  3X: 0.996134
  10X: 0.995013
  30X: 0.986197
  100X: 0.0785296
  300X: 0
  1000X: 0
```
Our minimum requirement is that the genome has 10X coverage over at least 95% of the genome. This sample has >70X average coverage, and more than 98% of the genome has >30X coverage, which is plenty for variant calling. Apart from the report itself (`DKC2-PP1.mapping-report.txt`), the program also generates a version of the reference genome `DKC2-PP1.masked-ref.fa.gz`, in which all positions with less than 5X coverage are masked out with Ns. This is useful for assessing gene absence/presence and to evaluate whether missing SNP-calls are caused by lack of coverage.

### Calling variants
We then use the script `call_variants.pl`, which uses a combination of the tools `mpileup` and `call` from Samtools to call variants on all bam-files stored in the sample folder
```bash
% scripts/call_variants.pl --reference reference/NC_000962.fna.gz --filter-regions reference/NC_000962.default_non-repetetive.bed example/DKC2-PP1
``` 
This creates three files: `DKC2-PP1.bcftools.raw.vcf.gz` (raw variant calls), `DKC2-PP1.bcftools.qfilt.vcf.gz` (quality filtered variants) and `DKC2-PP1.snps.vcf.gz` (only SNPs from non-repetetive parts of the genome)

### Creating a concatenated core-snp alignment
If we repeat the procedure on the other three sample runs, we should now have four sample folders, `DKC2-PP1`, `DKC2-PP15`, `DKC2-T101` and `DKC2-T91`. To make a core alignment of concatenated SNPs we use `vcf2fastas.pl`. This program scans folders for files with the extension `.snp.vcf.gz` and first compiles a list of 'high-quality' SNP-positions that lives up to certain minimum criteria involving, minimum read depth, minimum number of forward and reverse reads, and a minimum allele frequency. These criteria only need to be fulfilled in a single sample for a position to be included, and are adjustable through the options `--hq_mindepth (default: 10)`, `--hq_minfwd (default: 4)`, `--hq_minrev (default: 4)` and `--hq_minfreq (default: 0.85)`. Alternatively, it is possible to provide a bed-file of pre-defined regions, using the `--regions` option to be included in the alignment (note: disables the search for high quality positions). 

```bash
% scripts/vcf2fasta.pl --reference reference/NC_000962.fna.gz --snp-dist 12 --prefix example/core-snps example/DKC2-*
```

The option `--snp-dist` tells the program to filter out SNP-positions that are closer than 12 bp apart. It is usually a good idea to be suspicious of SNPs located very close to each other as, provided genomes are closely related, as these are often more likely to be a result of recombination or mis-mapped reads than stocahstichally determined point mutations.

After the scripts is done the following output files are generated:

Extension|Description
-----------------------------------|---------------------------------------------------
.alignment_positions.txt|An overview of which positions on the reference genome that were included in the final core alignment
.aln|The core alignment it self (in FASTA-format). Used as basis for phylogenetic analyses and calculation of pairwise distances. 
.excluded_positions.txt|Which positions of the genome that were excluded from the final core alignment during analysis for various reasons (e.g. monomorphic, ambiguous or too closely spaced SNP-positions)
.hq_postitions.bed|A bed-formatted file of all detected high-quality SNP-positions (not generated when using the `--regions` option)


### Generating a sample distance matrix
Finally, we generate a distance matrix with `aln2distmat.R`:
```
% scripts/aln2distmat.R example/core-snps.aln > example/core-snps.distmat.csv
```

Which should look like this:

Sample|DKC2-PP1|DKC2-PP15|DKC2-T101|DKC2-T91
:---------|:----------:|:----------:|:----------:|:----------:
DKC2-PP1|0|12|8|34
DKC2-PP15|12|0|12|42
DKC2-T101|8|12|0|36
DKC2-T91|34|42|36|0

# Code


## `scripts/pre_process_reads.sh`
```bash
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
```

## `scripts/map_reads.sh`
```bash
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
	# Locate existing files from the reference prefix
	# and list all found input files
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
```

## scripts/make_mapping_report.pl
```perl
#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

use Cwd;
use File::Basename;
use File::Find;
use File::stat;

use List::Util qw(sum max);

use Getopt::Long;
use Time::Piece;
use Time::Seconds;
use Digest::MD5 qw(md5_hex);

use Bio::SeqIO;

my $VERSION  = '0.4.11';

my %file_types = (
	'mapping'        => qr/^(.*)\.(sam|bam|cram)$/,
	'mapping-index'  => qr/^(.*)\.(bam\.bai|cram\.crai)$/,
	'depth'          => qr/^(.*)\.depth\.gz$/,
	'depth-index'    => qr/^(.*)\.depth\.gz\.(tbi)$/,
	'depth-report'   => qr/^(.*)\.depth-report\.txt$/,
	'mapping-stats'  => qr/^(.*)\.mapping-stats\.txt$/,
	'reference-mask' => qr/^(.*)\.masked-ref\.(fa|fna|fasta)\.gz$/,
	'unmapped-reads' => qr/^(.*)\.unmapped_R[12]\.(fq|fastq).gz$/
);

my $QUIET;
my $WARNINGS = 0;
my $INDENT   = " " x 3;

my $SEPARATOR = ": ";
my $WORKDIR = getcwd();

my $MAX_FIND_DEPTH      = 2;
my $MAX_LIBS            = 5;
my $MIN_READS           = 10_000;
my $MAX_MAPPED_ERR_RATE = 0.01;
my $MIN_MAPPED_FREQ     = 0.85;
my $MIN_MAPQ            = 30;
my $MIN_BAQ             = 20;
my $REF_MASK_MINDEPTH   = 5;
my $REF_MASK_UNKNOWN    = 'N';
my $REF_MASK_MISSING    = '-';

my @THRESH   = qw(
	0
	3
	10
	30
	100
	300
	1000
);

my ($BED_FILE, $REF_FILE, $INPUT_FILE);

&GetOptions(
	'quiet'         => \$QUIET,
	'regions|R=s'   => \$BED_FILE,
	'reference|T=s' => \$REF_FILE
);

my $dir = shift @ARGV;

(defined ($dir) && -d $dir) || err("Folder name ${dir} not valid: $!");

my %files       = ();
my %refseqs     = ();
my %refseqs_md5 = ();
my %regions     = ();
my %libraries   = ();

add_references($REF_FILE);
add_regions($BED_FILE);
scan_sample_folder($dir);

err("No valid data files found in folder $dir") unless $files{'files_found'};

#### MAPPING ####

my $mstat_files = $files{'mapping-stats_files'} || {};
my @mapped_libs = sort { $mstat_files->{$a}->{'mtime'} <=> $mstat_files->{$b}->{'mtime'} } keys %{ $mstat_files };
my @num_reads   = ();
my @map_freqs   = ();
my @err_rates   = ();

err("Too many (>${MAX_LIBS}) mapped libraries found in directory $dir $files{'sample'}") if @mapped_libs > $MAX_LIBS;
msg("Found mapping statistics for", scalar @mapped_libs, "librarires");

## Validate mapping statistics
foreach my $id(@mapped_libs) {
	my $mstat_file    = $mstat_files->{$id};
	my $mstats        = fetch_mapping_stats($mstat_file->{'file'});
	my $mapping_file  = $files{'mapping_files'}->{$id};
	my $mapping_index = $files{'mapping-index_files'}->{$id};

	my @mstat_keys = qw(
		sequences
		reads_mapped
		bases_mapped__cigar_
		error_rate
	); 

	# Validate number of sequences
	my ($seqs, $reads_mapped, $bases_mapped, $err_rate) = @{ $mstats }{@mstat_keys};
	
	msg("Checking library $id (mapping date: $mstat_file->{'mtime'})");
	msg($INDENT,"Filtered sequences:", $seqs, $seqs > $MIN_READS ? "OK" : "FAILED");

	push @num_reads, $seqs;

	# Validate mapping frequency
	my $map_freq = $reads_mapped / $seqs;
	msg($INDENT,"Mapping frequency:", 
		sprintf("%.2f%% (mapped bases: %i)", 100 * $map_freq, $bases_mapped),
		$map_freq > $MIN_MAPPED_FREQ ? "OK" : "FAILED");

	push @map_freqs, $map_freq;

	# Validate error rate
	msg($INDENT,"Error rate:", $err_rate, $err_rate < $MAX_MAPPED_ERR_RATE ? "OK" : "FAILED");

	push @err_rates, $err_rate;

	## Validate mapping files

	if ($mapping_file) {
		my ($fn,$mtime) = @{ $mapping_file }{qw/file mtime/};

		msg($INDENT, "Mapping file: $fn");

		if ($mapping_index) {
			if ($mapping_index->{'mtime'} >= $mtime) {
				msg($INDENT, "Index: $mapping_index->{'file'}", "OK");
			} else {
				msg($INDENT, "Index outdated ($mapping_index->{'mtime'}. Re-indexing $fn");
				system("samtools index $fn");
				msg($INDENT, "Done");
			}

		} else {
				msg($INDENT, "Indexing $fn");
				system("samtools index $fn");
				msg($INDENT, "Done");
		}

		# Validate mapping-file header
		my $bam_header = fetch_bam_header($mapping_file->{'file'});

		msg($INDENT, "Validating header");

		# Validate mapping-file reference sequences
		my @SQ = @{ $bam_header->{'@SQ'} } || ();

		foreach my $seq(@{ $bam_header->{'@SQ'} }) {
			my ($seqid, $length, $md5) = @{ $seq }{qw/SN LN M5/};

			if (defined($seqid) && exists $refseqs{$seqid}) {
				err("Reference sequence $seqid in BAM-header has a different length than in $REF_FILE")
					unless $refseqs{$seqid}->{'len'} == $length;
				err("Sequence checksums (MD5) for $seqid do not match between BAM-header and $REF_FILE")
					if (defined($md5) && $refseqs{$seqid}->{'md5'} ne $md5);
				msg($INDENT, "Reference sequence $seqid OK")
			} elsif (defined($md5) && exists $refseqs{$md5}) {
				wrn("Sequence identifier ($seqid) not found in $fn, but but an identical reference sequence has been found in $REF_FILE ($refseqs{$md5}->{'id'})");
			} else {
				wrn("Sequence with ID $seqid present in $fn not found in $REF_FILE");
			}
	 	}

	 	# Validate read groups
	 	my @RG = exists $bam_header->{'@RG'} ? @{ $bam_header->{'@RG'} } : ();

	 	if (@RG == 1) { 
	 		my ($rg_id, $rg_sample) = @{$RG[0]}{qw/ID SM/};

	 		if ($rg_sample eq $files{'sample'}) {
	 			msg ($INDENT, "Sample ID: $rg_sample","OK");
	 		} else {
	 			wrn("Read group sample name ($rg_sample) does not match sample folder name ($files{'sample'})");
	 		}

	 		if ($rg_id eq $id) {
	 			msg($INDENT, "Read group: $rg_id", "OK");
	 		} else {
	 			wrn("Read group in $fn contains unkown identifier ($id)");
	 		}
	 	} elsif (@RG > 1) {
	 		wrn("Mapping file $fn has more than one read group");
	 	} else {
	 		wrn("No read group information for mapping file $fn");
		}
	}
}

#### COVERAGE ####

# Check depth file
my $depth_file = $files{'depth_files'}->{ $files{'sample'} };

if (defined ($depth_file)) {
	my ($fn,$mtime,$id) = @{$depth_file}{qw/file mtime id/};

	msg("Validating depth file $fn (scan date: $mtime)");
	my $depth_index = $files{'depth-index_files'}->{$id};

	if (defined ($depth_index) && $depth_index->{'mtime'} >= $mtime) {
		msg($INDENT,"Index $depth_index->{'file'}","OK");
	} else {
		msg($INDENT, "Indexing $fn");
		system("tabix $fn -s 1 -b 2 -e 2");
		msg($INDENT, "Done");
	}

	my $TAB;

	# Verify sequence names in depth index
	open ($TAB, "tabix -l $fn|");
	chomp(my @seqids = <$TAB>);
	close $TAB;

	foreach my $seqid(@seqids) {
		if (exists $refseqs{$seqid}) {
			msg($INDENT, "$seqid", "OK");
		} else {
			wrn("$seqid not present in $REF_FILE");
		}
	}

	# Verify number of libraries
	open ($TAB, "bgzip -cd $fn|");
	chomp(my $firstline = <$TAB>);
	
	my @F = split /\t/, $firstline;
} else {
	msg("Scanning mapping files for depth");

	my $fn = $files{'folder'} . "/$files{'sample'}.depth.gz";
	my @map_files = map { $_->{'file'} } @{ $files{'mapping_files'} }{@mapped_libs};

	system("samtools depth --reference $REF_FILE -a -Q $MIN_MAPQ -q $MIN_BAQ @map_files | bgzip > $fn");
	system("tabix $fn -s 1 -b 2 -e 2");

	msg("Done");

	$depth_file->{'file'} = $fn;
}

### MASK ###

# Check reference mask
my $refmask_file = $files{'reference-mask_files'}->{ $files{'sample'} };

if (defined ($refmask_file)) {
	msg("Validating reference mask $refmask_file->{'file'}");

	open (my $FH, "bgzip -cd $refmask_file->{'file'}|");

	my $mask_fh = Bio::SeqIO->new( -fh => $FH, -format => 'Fasta');

	while (my $seq = $mask_fh->next_seq) {
		my $seqid  = $seq->id;
		my $length = $seq->length;

		err("$seqid not found") unless (exists $refseqs{$seqid});
		msg($INDENT, "$seqid OK");
		err("length does not match") unless ($refseqs{$seqid}->{'len'} == $length);
		msg($INDENT, "Length: $length bp OK");

		my $nt = $seq->seq =~ tr/ACTG/ACTG/;

		msg($INDENT, "Coverage OK (" . sprintf("%.2f%%", (100 * $nt) / $length) . ")");
	}

	close $FH;
} else {
	msg("Generating reference mask");

	my %masked = ();

	foreach my $id(keys %refseqs) {
		$masked{$id} = $REF_MASK_MISSING x $refseqs{$id}->{'len'};
	}

	open (my $FH, "bgzip -cd $depth_file->{'file'} |");
	while (defined (my $line = <$FH>)) {

		chomp $line;

		my ($seqid,$pos, @depth) = split /\t/, $line;
		my $cov = sum(@depth);

		next unless $cov > 0;

		my $nt = $cov < $REF_MASK_MINDEPTH ? $REF_MASK_UNKNOWN : substr($refseqs{$seqid}->{'seq'}, $pos - 1, 1);

		substr($masked{$seqid}, $pos - 1, 1, $nt);
	}

	close $FH;

	my $refmask_fn = "$files{'folder'}/$files{'sample'}.masked-ref.fa.gz";
	open ($FH, "| bgzip -c > $refmask_fn");

	msg($INDENT, "Writing to $refmask_fn");
	my $fasta_fh = Bio::SeqIO->new( -format => 'fasta', -fh => $FH);

	for my $id(sort keys %masked) {
		$fasta_fh->write_seq( Bio::Seq->new(-id => $id, -seq => $masked{$id}));
		msg($INDENT,"$id");
	}

	close $FH;

}

### STATISTICS ###

my %histogram = ();
my $sample    = $files{'sample'};
my $sum       = 0;
my $num_bams  = 0;
my ($min,$max);

my %counts    = ();
do { $counts{$_} = 0 } for @THRESH;

msg("Generating coverage report");

my $depth_fn = $depth_file->{'file'};

open (my $TAB, ((defined ($BED_FILE) && -r $BED_FILE) ? "tabix -R $BED_FILE $depth_fn |" : "bgzip -cd $depth_fn |"));

while (defined (my $line = <$TAB>)) {
	chomp $line;
	
	my @F   = split /\t/,$line;
	
	splice(@F,0,2);
	
	my $cov = sum(@F);

	if ($. == 1) {
		$num_bams  = scalar @F;
		$min = $cov;
		$max = $cov;
	}

	$min  = $cov < $min ? $cov : $min;
	$max  = $cov > $max ? $cov : $max;
	$sum += $cov;
	
	$histogram{$cov}++;

	do { $counts{$_}++ if $cov >= $_ } for @THRESH;
}

close $TAB;

my %result = ();
my $base = $counts{ shift @THRESH };

err("No bases mapped to reference!") unless $base;

my $sigma = 0;
my $mean  = $sum / $base;

while (my ($cov,$val) = each %histogram) {
	$sigma += $val * ($cov - $mean)**2;
}

my $stdev = sqrt($sigma / $base);
my $outsidelim = 0;
my $minlim = max(0, $mean - (2 *$stdev));
my $maxlim = $mean + (2* $stdev);

while (my ($cov,$val) = each %histogram) {
	$outsidelim += ($cov * $val) if ($cov < $minlim || $cov > $maxlim);
}

my $rating = 0;

do { $rating++ if $_ / $base >= 0.95} for @counts{@THRESH};

my @info_keys = qw(
	version
	sample
	number_of_libraries
	library_names
	library_number_of_reads
	library_mapping_frequencies
	library_error_rates
	bases_mapped
	positions_covered
	average_coverage
	stdev_coverage
	coverage_outside_2sd
	minimum_coverage
	maximum_coverage
	coverage_rating
	coverage_thresholds
);

@result{@info_keys} = (
		$VERSION,
		$sample,
		$num_bams,
		join(",", @mapped_libs),
		join(",", @num_reads),
		join(",", @map_freqs),
		join(",", @err_rates),
		$sum,
		$base,
		$mean,
		$stdev,
		$outsidelim / $sum,
		$min,
		$max,
		$rating,
		''
);

foreach my $key(@info_keys) {
	print join( $SEPARATOR, $key, $result{$key} ) . "\n";
}

foreach my $thr(@THRESH) {
	my $key = $thr . "X";
	print "  " . join( $SEPARATOR, $key, sprintf("%g", $counts{$thr} / $base )) . "\n"; 
}

msg("Done. There were $WARNINGS warnings");


##### SUBROUTINES #####

sub scan_sample_folder {
	my ($dir) = @_;

	return unless -d $dir;

	$dir = $dir eq '.' ? getcwd() : $dir;

	my $sample = dir_to_id($dir);

	msg("Scanning sample folder $sample (path: $dir)");
	
	$files{'sample'} = $sample;
	$files{'folder'} = $dir;
	$files{'root_depth'} = $dir =~ tr[/][];

	find( { wanted => \&add_files, preprocess => \&check_dir_depth }, $dir );

	msg("$files{'files_found'} data files found");
}

sub check_dir_depth {
	my $depth = ($File::Find::dir =~ tr[/][]) - $files{'root_depth'};

	return @_ if $depth < $MAX_FIND_DEPTH;
	return grep { not -d } @_ if $depth == $MAX_FIND_DEPTH;
	return;
}

sub add_files {
	return unless -r -f $_;

	my $fn    = $_;
	my $real  = $File::Find::name;

	while (my ($tag,$regex) = each %file_types) {
		if ($fn =~ m/$regex/) {
			my $key   = $tag . "_files";
			my $stat  = stat($fn);
			my $mtime = scalar localtime $stat->mtime;
			my $id    = basename($1);
			my $type  = $2;
			my %info = ( file => $real, mtime => $mtime, id => $id, type => $type );
			
			msg($INDENT,"Found $tag file $fn");

			if (exists $files{$key}->{$id} ) {
				if (ref($files{$key}->{$id}) eq 'ARRAY') {
					push @{ $files{$key}->{$id} }, \%info;
				} else {
					$files{$key}->{$id} = [ $files{$key}->{$id}, \%info ];
				}
			} else {
				$files{$key}->{$id} = \%info;
			}

			$files{'files_found'}++;
		}
	}
}

sub add_references {
	my ($file) = @_;
	return () unless (defined $file && -r $file);

	my @seqs   = ();

	open ( my($FH), $file =~ m/\.gz$/ ? "bgzip -cd $file |" : "$file");

	my $seqin = Bio::SeqIO->new( -fh => $FH, -format => 'Fasta');
	msg("Loading reference sequence from $REF_FILE");

	while (my $seq_o = $seqin->next_seq) {
		push @seqs, { 
			id  => $seq_o->id,
			seq => uc($seq_o->seq),
			len => $seq_o->length,
			md5 => md5_hex($seq_o->seq)
		}
	}

	close $FH;

	my $reflen = 0;
	my $num_refseqs = scalar @seqs;
	do { $reflen += $_->{'len'}; } for @seqs;

	msg("$num_refseqs sequences added (length: $reflen)");

	%refseqs     = map { $_->{'id'} => $_  } @seqs;
	%refseqs_md5 = map { $_->{'md5'} => $_ } @seqs;
}

sub add_regions {
	my ($file) = @_;

	return unless (defined $file && -r $file);

	$files{'bed_file'} = $file;

	msg("Loading regions from $file");

	open (my $BED, "<", $file);

	my $totlen  = 0;
	my $numregs = 0;

	while (defined (my $line = <$BED>)) {
		chomp $line;

		my ($chrom, $beg, $end) = split /\t/, $line, 3;

		$totlen += ($end - $beg);
		$numregs++
	}

	msg("$numregs regions covering $totlen bases added");
}

sub dir_to_id {
	my($dir) = @_;
	$dir =~ s{/$}{}; # remove trailing slash first
	
	my @dir = File::Spec->splitpath($dir);

	return $dir[-1];  
}

sub fetch_bam_header {
	my ($bam_file) = @_;

	my %header = ();

	open (my $BAM, "samtools view -H $bam_file|");

	while (defined (my $line = <$BAM>)) {
		chomp $line;

		my @F = split "\t", $line;

		my $key = shift @F;

		last unless $key =~ m/^@[A-Z]{2}$/;

		push @{ $header{$key} }, { map { ($_ =~ m/^([A-Z0-9]{2}):(.*)$/) } @F };
	}
	
	return \%header;
}

sub fetch_mapping_stats {
	my ($stats_file) = @_;

	open (my $STAT, "<", $stats_file);

	my %stats = ();

	while (defined (my $line = <$STAT>)) {
		next if $line =~ m/^#/;
		chomp $line;

		if ($line =~ m/^SN\t([^\t]+):\t([^\t]+)/) {
			my ($key,$val) = ($1,$2);

			$key =~ tr/[ ()]/_/;

			$stats{$key} = $val;
		}
	}

	return \%stats;
}

sub msg {
	return if $QUIET;
	my $t = localtime;
	my $string = "[" . $t->hms . "] @_\n";

	print STDERR "$string";
}

sub wrn {
	$WARNINGS++;
	return if $QUIET;
	msg("WARNING:", @_);
}

sub err {
	$QUIET=0;
	msg("ERROR:", @_);
	exit(2);
}
```

## `scripts/call_variants.pl`

```perl
#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

use Cwd;
use Getopt::Long;

use File::Basename;
use File::Spec;
use File::Slurp;

use Bio::SeqIO;

my ($PREFIX,$BED);

my $MINMAPQ   = 30; # Minimum mapping quality (MapQ)
my $MINBAQ    = 13; # Minumum base quality (BAQ)
my $MINIREADS = 5;  # Minimum number of gapped reads for indel-candidates
my $MINQUAL   = 20; # Minimum VCF-qual
my $MINCOV    = 2;  # Minimum number of both forward and reverse reads needed to retain a raw snp call

my $CWD       = getcwd();
my $REF       = $ENV{'SANDBOX_REF_PATH'};
my $REF_FILT  = $ENV{'SANDBOX_FILT_PATH'};

my @mpileup_tags = qw(
	DP
	AD
	ADF
	ADR
	SP
);

&GetOptions(
	'prefix=s'         => \$PREFIX,
	'reference=s'	     => \$REF,
	'filter-regions=s' => \$REF_FILT,
	'pileup-regions=s' => \$BED

);

my $dir = shift @ARGV;

(-d $dir) or die ("Cannot read from $dir, not a valid directory");

$dir = File::Spec->rel2abs($dir);

# Check reference sequence
(-r $REF) or die("Cannot read from reference file $REF");

warn "Using reference $REF\n";

$REF      = File::Spec->rel2abs($REF);
$REF_FILT = File::Spec->rel2abs($REF_FILT) if $REF_FILT;

#open (my $REFIN, $REF =~ m/\.gz$/ ? "bgzip -cd $REF |" : "$REF");

my $seqin      = Bio::SeqIO->new( -fh => $REFIN, -format => 'fasta');
my $reflen     = 0;
my @refseq_ids = ();

while (my $seq = $seqin->next_seq) {
	my $seqid = $seq->id;

	$reflen += $seq->length;

	push @refseq_ids, $seqid;
}

my $num_refseqs = scalar @refseq_ids;

close $REFIN;

warn "\n$num_refseqs reference sequence(s) found (length: $reflen bp)\n";

if ($BED) {
	$BED = File::Spec->rel2abs($BED);

	open (my $BEDIN, "<", $BED);

	my $npos = 0;
	warn "\nLoading pre-defined regions from $BED\n";
	while (defined (my $line = <$BEDIN>)) {
		chomp $line;

		my @F = split /\t/,$line;
		
		die("Wrongly formatted BED-file") unless (@F >= 3);

		$npos += ($F[2] - $F[1]);
	}

	my $npos_per = sprintf("%.4f%%", (100 * $npos)/ $reflen);
	warn "\nVariants will be called at $npos ($npos_per) positions\n";
}

(-d $dir) or die("Input directory $dir does not exist");

$PREFIX = dir_to_id($dir) unless $PREFIX;

warn "\nCalling variants in directory $dir, using prefix $PREFIX\n";

opendir(my $DIR, $dir);

my @commands = ();
my @bams = grep { /\.(?:bam|cram)$/ } readdir $DIR;

warn "\nFound following mapping files: ", join(", ",@bams),"\n\n";

my $log_file = "${PREFIX}.variant_caller.log";

# Generate variant caller pipe
my $mpileup_cmd = join(" ", 
	"samtools mpileup",
	"-m $MINIREADS",
	($BED ? "-l $BED" : ()),
	"-q $MINMAPQ",
	"-Q $MINBAQ",
	"-uf $REF -t " . join (",",@mpileup_tags), 
	@bams
);

my $call_cmd = join (" ", 
	"bcftools call",
	($BED ? "-mA" : "-mvA"),
	"--ploidy 1"
);

my $norm_cmd = "bcftools norm -f $REF -O z -o ${PREFIX}.bcftools.raw.vcf.gz";

push @commands, "(" . join (" | ", $mpileup_cmd, $call_cmd, $norm_cmd) . ")";
push @commands, "bcftools filter -o ${PREFIX}.bcftools.qfilt.vcf.gz -O z -i \'QUAL >= $MINQUAL & (DP4[2] >= $MINCOV & DP4[3] >= $MINCOV)\' ${PREFIX}.bcftools.raw.vcf.gz";
push @commands, "(bcftools view -v snps " . ($REF_FILT ? "-T $REF_FILT" : "") . " -O z -o ${PREFIX}.snps.vcf.gz ${PREFIX}.bcftools.qfilt.vcf.gz)";

chdir("$dir");

foreach my $cmd(@commands) {
	warn "Running: $cmd\n";
	
	append_file( $log_file, "\n### $cmd\n\n" );

	$cmd .= " 2>> $log_file";

	system("$cmd")== 0 or die("Error running command");
}

chdir("$CWD");

sub dir_to_id {
	my($dir) = @_;
	$dir =~ s{/$}{}; # remove trailing slash first
	
	my @dir = File::Spec->splitpath($dir);

	return $dir[-1];  
}
```

## scripts/vcf2fasta.pl

```perl
#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

use File::Basename;
use File::Spec;

use Bio::SeqIO;
use Bio::AlignIO;
use Bio::LocatableSeq;

use Getopt::Long;
use Vcf;

use List::Util qw(sum);

my %hq_pos        = ();
my %core_pos      = ();
my %refseq_seq    = ();
my %vcf_files     = ();
my @refseq_ids    = ();
my @sample_ids    = ();

my $MINDEPTH     = 10;
my $MINFWD       = 4;
my $MINREV       = 4;
my $MINFREQ      = 0.85;
my $MIN_SNP_DIST = 2;
my $NOFILTER;
my $FH;

my $INBED    = '';
my $INPREFIX = 'snps';
my $PREFIX   = 'core';	
my $CHROM    = '';#$ENV{'NCBI_H37RV'};
my $REF      = $ENV{'SANDBOX_REF_PATH'};

&GetOptions(
	'hq_mindepth=i'  => \$MINDEPTH,
	'hq_fwd=i'       => \$MINFWD,
	'hq_rev=i'       => \$MINREV,
	'hq_freq=f'      => \$MINFREQ,
	'regions|R=s'    => \$INBED,
	'disable-filter' => \$NOFILTER,
	'snp-dist=i'     => \$MIN_SNP_DIST,
	'inprefix=s'     => \$INPREFIX,
	'prefix=s'       => \$PREFIX,
	'reference|r=s'  => \$REF
);

open ($FH, $REF =~ m/\.gz$/ ? "bgzip -cd $REF |" : "$REF");

my $seqin = Bio::SeqIO->new( -fh => $FH, -format => 'Fasta');
my $reflen = 0;

warn "Loading reference sequence from $REF\n";

while (my $seq = $seqin->next_seq) {
	my $id = $seq->id;

	push @refseq_ids,$id;

	$refseq_seq{ $id } = uc($seq->seq);
	$reflen 	  += $seq->length;
}

my $num_refseqs = scalar @refseq_ids;

close $FH;

warn "$num_refseqs sequences loaded (length: $reflen)\n";
warn "Sniffing sample folders...\n";

my $num_samples = 0;
my $num_vcfs    = 0;

foreach my $dir(@ARGV) {
	my $id = &dir_to_id($dir);
	my $vcf_file = "$dir/$id.$INPREFIX.vcf.gz";

	if ( -e "$vcf_file" ) {
		print STDERR "\rFound $vcf_file" . " " x (80 - length($vcf_file));
		my $vcf = Vcf->new( file => "$vcf_file" );

		$vcf->parse_header();

		my @samples = $vcf->get_samples();

		$num_vcfs++;
		$num_samples += scalar @samples;
		$vcf->close();

		$vcf_files{ $vcf_file } = \@samples;
	}
}

print STDERR "\n$num_samples samples found.\n";

if ($INBED) {
	open ($FH, "<", $INBED);

	while (defined (my $line = <$FH>)) {
		chomp $line;
		my @F = split /\t/, $line;
		$hq_pos{$F[0]}->{$F[1] + 1}->{$F[3]} = 1
	}

	close $FH;
} else {
	warn "\nScanning for high quality snp positions...\n";

	foreach my $vcf_file(sort keys %vcf_files) {
		print STDERR "\rLooking at $vcf_file" . " " x (80 - length($vcf_file));

		my $vcf = Vcf->new( file => "$vcf_file" );

		$vcf->parse_header();

		my @samples = @{ $vcf_files{$vcf_file}};

		while (my $v  = $vcf->next_data_hash()) {
			my $pos   = $v->{'POS'};
			my $info  = $v->{'INFO'};
			my $seqid = $v->{'CHROM'};

			next if exists $info->{'INDEL'};

			foreach my $sample(@samples) {
				die ("_REFERENCE is an illegal sample name!") if $sample eq '_REFERNCE';

				my $gtype = $v->{'gtypes'}->{$sample};

				my $dp  = $gtype->{'DP'};
				my @ad  = split/,/, $gtype->{'AD'};
				my @adf = split/,/, $gtype->{'ADF'};
				my @adr = split/,/, $gtype->{'ADR'};

				next unless
					sum(@ad)  >= $MINDEPTH &&
					sum(@adf) >= $MINFWD &&
					sum(@adr) >= $MINREV;

				my $freq = $ad[1] / sum(@ad);
				my $ref  = $v->{'REF'};
				my $alt  = $v->{'ALT'}->[0];

				$hq_pos{$seqid}->{ $pos }->{ "${ref}->${alt}" }++; 
			}
		}

		$vcf->close();
	}
}

my %remove_snps = map { $_ => {} } @refseq_ids;
my %positions = ();
my %templates = ();
my $num_pos = 0;

foreach my $seqid(@refseq_ids) {
	my $pos_r = $hq_pos{$seqid} || {};

	$num_pos += scalar keys %{ $pos_r };
	$positions{$seqid}= [ sort { $a <=> $b } keys %{$pos_r} ];
}

print STDERR "\nA total of $num_pos genome positions will be investigated.\n";

open ($FH, '>', "$PREFIX.hq_postitions.bed");

foreach my $seqid(@refseq_ids) {
	foreach my $pos(@{ $positions{$seqid} }) {
		print $FH join("\t", $seqid, $pos - 1, $pos, each %{ $hq_pos{$seqid}->{$pos} }) . "\n";
	}
}

close $FH;

# Scan for closely spaced snps and filter them out
if ($MIN_SNP_DIST && ! $NOFILTER) {
	foreach my $seqid(@refseq_ids) {
		my @positions = @{ $positions{$seqid} };
		
		next unless @positions >= 2;

		my $i  = $#positions + 1;

		while ( --$i >= 0 ) {
			my $pos   = $positions[ $i     ];
			my $left  = $positions[ $i - 1 ]  if $i > 0;
			my $right = $positions[ $i + 1 ]  if $i < $#positions;

			my $dist = $MIN_SNP_DIST;

			if (defined($left)) {
				$dist = $pos - $left;
			} 

			if (defined($right)) {
				my $dist_r = $right - $pos;
				$dist = $dist_r < $dist ? $dist_r : $dist;
			}

			$remove_snps{$seqid}->{$pos} = { 'reason' => 'TOO_CLOSE', 'stats' => { '-' => $num_samples } }
				if $dist < $MIN_SNP_DIST;
		}
	}
}

# Construct a universal reference template
foreach my $seqid(@refseq_ids) {
	my @positions = grep { ! exists $remove_snps{$seqid}->{$_} }@{ $positions{$seqid} };

	foreach my $pos(@positions) {
		$templates{'_REFERENCE'}->{ $seqid } .= substr($refseq_seq{$seqid}, $pos - 1, 1); 
	}

	$positions{$seqid} = \@positions;
}

# Collect reference templates from all samples
warn "Collecting sample masks...\n";

foreach my $vcf_file(sort keys %vcf_files) {
	my $folder = dirname($vcf_file);

 	foreach my $sample(@{ $vcf_files{$vcf_file} }) {
 		if (-e "$folder/$sample.masked-ref.fa.gz") {
	 		open (my $FH, "bgzip -cd $folder/$sample.masked-ref.fa.gz|");

	 		my $mask_fh = Bio::SeqIO->new( -fh => $FH, -format => 'Fasta');
	 		my $nts = 0;

			while (my $seq = $mask_fh->next_seq) {
				my $template = '';
				my $mask     = $seq->seq;
				my $seqid    = $seq->id;

				$nts += $mask =~ tr/ACTG/ACTG/;

				my @positions = @{ $positions{$seqid} };

				do { $template .= substr($mask, $_ - 1, 1) } for @positions;

				$templates{ $sample }->{$seqid} = $template;	
			}

			close $FH;

			my $covered = sprintf("%.4f", 100 * $nts / $reflen);

			warn "\nWARNING: sample $sample is poorly covered (${covered}%)\n" if ($covered < 95);

			print STDERR "\r$sample (${covered}% covered. " . --$num_samples . " left)   ";
		} else {
			print STDERR "\r$sample (no mask found. "       . --$num_samples . " left)  \n";
		}
	}
}

warn "\nCompiling alignment...\n";

my %pos_stats  = ();
my %pos_freqs  = ();
my %aln_seqs   = ();
my %aln_offset = ();

foreach my $vcf_file(sort keys %vcf_files) {

	my $vcf = Vcf->new( file => $vcf_file);
	
	$vcf->parse_header();

	my @samples = @{ $vcf_files{$vcf_file} };
	my %gts = ();

	while (my $v = $vcf->next_data_hash()) {
		next if exists $v->{'INFO'}->{'INDEL'};

		my $pos   = $v->{'POS'};
		my $chrom = $v->{'CHROM'};

		next unless (exists $hq_pos{$chrom}->{$pos} && ! exists $remove_snps{$chrom}->{$pos});

		foreach my $sample(@samples) {
			my $gt = $v->{'gtypes'}->{$sample};
			my @ad = split/,/, $gt->{'AD'};

			next unless (@ad > 1 && sum(@ad));

			my $freq = $ad[1] / sum(@ad);

			push @{ $pos_freqs{$chrom}->{$pos}}, sprintf("%.4f", $freq);

			$gts{$chrom}->{$pos}->{$sample} = $v->{'ALT'}->[ $gt->{'GT'} - 1 ] ;
		}
	}
	
	$vcf->close();

	foreach my $sample(@samples) {
		$aln_seqs{$sample} = '';

		foreach my $seqid(@refseq_ids) {
			my $template = exists $templates{$sample}->{$seqid} ? 
				$templates{$sample}->{$seqid} :
				$templates{'_REFERENCE'}->{$seqid};

			my $seq = '-' x length($template);
			
			my $i = -1;
			my @positions = @{ $positions{$seqid} };

			die("Length mismatch") unless @positions == length($seq);

			while ( ++$i < length($seq) ) {

				my $pos = $positions[$i];
				
				my $nt = 'N';
				my $tmp_nt = substr($template, $i, 1);

				if (($tmp_nt ne 'N' && $tmp_nt ne '-') && 
					exists $gts{$seqid}->{$pos} && 
					exists $gts{$seqid}->{$pos}->{$sample}) 
				{
					$nt = $gts{$seqid}->{$pos}->{$sample};
				} else {
					$nt = $tmp_nt;
				}

				$pos_stats{$seqid}->[ $i ]->{$nt}++;

				substr($seq, $i, 1, $nt);
			}

			$aln_offset{$seqid} = length($aln_seqs{$sample});
			$aln_seqs{$sample} .= $seq;
		}

		print STDERR "\r$sample";
		push @sample_ids, $sample;
	}
}

$num_samples = scalar keys %aln_seqs;

my $unknown   = 0;
my $mono      = 0;

unless ($NOFILTER) {
	warn "\nScanning alignment for junk...\n";

	foreach my $seqid(@refseq_ids) {
		my @positions = @{ $positions{$seqid} };
		my @pos_stats = @{ $pos_stats{$seqid} };

		my $i = @positions;

		while ( --$i >= 0) {
			my $pos      = $positions[$i];
			my $stat     = $pos_stats[$i];

			my @nts      = grep { /[ACTG]/ } keys %{ $stat };
			my $called   = sum(@{ $stat }{@nts}) || 0;
			my $call_frq = $called / $num_samples;
			my $reason   = '';

			if (@nts == 1 && $call_frq >= 0.95) {
				$reason = 'MONOMORPHIC';
				$mono++;
			} elsif ($call_frq < 0.95) {
				my $gaps =   exists $stat->{'-'} ? $stat->{'-'} : 0;
				my $ambigs = exists $stat->{'N'} ? $stat->{'N'} : 0;

				$reason = $gaps >= $ambigs ? 'GAP' : 'AMBIGUOUS';
				$unknown++;	
			}

			if ($reason) {
				$remove_snps{$seqid}->{$pos} = { 'reason' => $reason, 'stats' => $stat };
				splice(@positions, $i, 1);
				splice(@pos_stats, $i, 1);

				foreach my $sample(keys %aln_seqs) {
					my $removed = substr($aln_seqs{$sample},$aln_offset{$seqid} + $i,1,'');
				}
			}

		}

		@{ $positions{$seqid} } = @positions;
		@{ $pos_stats{$seqid} } = @pos_stats;
	}

	warn "Removed $unknown ambiguous and $mono monomorphic positions...\n";

	open (my $EXCL, ">", "$PREFIX.excluded_positions.txt");

	print $EXCL join("\t", qw/CHROMOSOME POS VARIANT A T C G N - REASON/), "\n";

	foreach my $seqid(@refseq_ids) {
		foreach my $pos(sort { $a <=> $b } keys %{ $remove_snps{$seqid} }) {
			my ($var) = keys %{ $hq_pos{$seqid}->{$pos} };
			my $h_ref = $remove_snps{$seqid}->{$pos};
		
			my ($stat,$reason) = @{ $h_ref }{qw/stats reason/};
			my @fields = map { defined $_ ? $_ : 0 } @{ $stat }{qw/A T C G N -/};

			print $EXCL join("\t", $seqid, $pos, $var, @fields, $reason), "\n";

		}
	}
}

open (my $INCL, ">", "$PREFIX.alignment_positions.txt");
print $INCL join("\t", qw/CHROMOSOME POS VARIANT ALN_POS A T C G N - /), "\n";

foreach my $seqid(@refseq_ids) {
	my @positions = @{ $positions{$seqid} };
	my @pos_stats = @{ $pos_stats{$seqid} };

	my $i = -1;

	while ( ++$i < @positions) {
		my $pos      = $positions[$i];
		my $stat     = $pos_stats[$i];
		my ($var) = keys %{ $hq_pos{$seqid}->{$pos} };
		my @fields = map { defined $_ ? $_ : 0 } @{ $stat }{qw/A T C G N -/};

		print $INCL join("\t", $seqid, $pos, $var, $i + 1, @fields ), "\n";
	}
}

my $aln_fn  = "$PREFIX.aln";
my $aln_obj = Bio::SimpleAlign->new( -id => $aln_fn );

foreach my $sample(@sample_ids) {
	$aln_obj->add_seq( Bio::LocatableSeq->new( -id => $sample, -seq => $aln_seqs{$sample}, -start => 1 ) );
}

my $aln_out = Bio::AlignIO->new( -file => ">$aln_fn", -format => "fasta" );

warn "Final alignment is " . $aln_obj->length . " bp\n";

$aln_obj->set_displayname_flat();
$aln_out->write_aln($aln_obj);

sub dir_to_id {
	my($dir) = @_;
	$dir =~ s{/$}{}; # remove trailing slash first
	
	my @dir = File::Spec->splitpath($dir);

	return $dir[-1];  
}

sub num2range {
  local $_ = join ',' => sort { $a <=> $b } @_;
  s/(?<!\d)(\d+)(?:,((??{$++1})))+(?!\d)/$1-$+/g;
  return $_;
}
```

*Last modified: 2019-02-14*
