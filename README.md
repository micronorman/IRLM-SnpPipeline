# IRLM-SnpPipeline

NGS pipeline for calling variants (SNPs) of *Mycobacterium tuberculosis* at IRLM Copenhagen, Denmark

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

*Last modified: 2019-02-19*
