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