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
	'reference=s'	   => \$REF,
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

open (my ($REFIN), $REF =~ m/\.gz$/ ? "bgzip -cd $REF |" : "$REF");

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
