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
