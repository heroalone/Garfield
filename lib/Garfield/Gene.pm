package Garfield::Gene;

use strict;
use warnings;
use File::Basename;
use FindBin;
use Pod::Usage;
use Cwd qw(abs_path);
use Garfield::CommonArgument;
use Fcntl qw(:flock LOCK_EX LOCK_UN);
use Parallel::ForkManager;
use File::Spec;

my $max_retries = 10;
my $maf_filter = 0.02;
my ($out_tped_fh, $out_dnf_fh);

sub new {
	my $class = shift;
	return bless {}, $class;
}

sub check_para {
	my ($class, $opts) = @_;
	my $def = 0;
	my $num = 0;
	foreach (values %$opts) {
		if (defined $_) {
			if (!/ARRAY/ || @{$_}) {
				$def++;
			}
		}
		$num++;
	}
	if ($def == 0) {
		return 0;
	} else {
		return 1;
	}
}

sub check_para_sub {
	my ($class, $opts_sub, $opts) = @_;
	if ($opts->{help}) {
		pod2usage(-exitval => 0, -verbose => 2, -input => "$FindBin::Bin/doc/Garfield_Gene.txt");
		exit 0;
	}

	if (!check_para($class, $opts_sub)) {
		print "Required parameters missing!\n";
		pod2usage(-exitval => 1, -verbose => 2, -input => "$FindBin::Bin/doc/Garfield_Gene.txt");
		exit 0;
	}

	my $exit_code = 0;

	if (!$opts_sub->{"bed"}) {
		print "Please provide --bed file!\n";
		$exit_code++;
	} else {
		if (!-e $opts_sub->{"bed"}) {
			print "bed file: $opts_sub->{bed} does not exist. Please check!\n";
			$exit_code++;
		}
		$opts_sub->{"bed"} = abs_path $opts_sub->{"bed"};
	}


	#### Common arguments
	my $cm_arg = Garfield::CommonArgument->new();
	my $exit_num_return = $cm_arg->common_argument($opts_sub, $exit_code);
	$exit_code += $exit_num_return;

	if ($exit_code > 0) {
		pod2usage(-exitval => 1, -verbose => 2, -input => "$FindBin::Bin/doc/Garfield_Gene.txt");
		exit 0;
	} else {
		# return "TRUE";
		print "Gene bed file is: $opts_sub->{bed}\n";
		print "With the extension of flanking of each gene: $opts_sub->{extension} bp\n";
		print "\#threads will be used: $opts_sub->{'threads'}\n";
		RUN_Gene($class, $opts_sub);
	}
}

sub RUN_Gene {
	my ($class, $opts_sub) = @_;
	my $genotype_file = $opts_sub->{genotype_file};
	my $trait_file = $opts_sub->{trait_file};
	my $max_threads = $opts_sub->{threads};
	my $bed_file = $opts_sub->{bed};
	my $extension = $opts_sub->{extension};
	my $prefix = $opts_sub->{prefix};
	my $temporary = $opts_sub->{temporary};
	my $output_dir = $opts_sub->{outdir};
	my $keep_negative = $opts_sub->{keep_negative};
	
	open($out_tped_fh, '>', "$output_dir/Garfield.Geno.$prefix.tped") or die "Cannot open $output_dir/Garfield.Geno.$prefix.tped for writing: $!";
	open($out_dnf_fh, '>', "$output_dir/Garfield.bestDNF.$prefix.txt") or die "Cannot open $output_dir/Garfield.bestDNF.$prefix.txt for writing: $!";
	

	# my (%TPED, %bestDNF);
	my $pm = Parallel::ForkManager->new($max_threads);
	
	open(my $bed_fh, '<', $bed_file) || die "Cannot open $bed_file: $!";
	while (my $line = <$bed_fh>) {
		chomp $line;
		$line =~ s/\r//g;
		next if $line =~ /^#/ || $line =~ /^\s*$/;

		my @i = split("\t", $line);
		if (@i < 4) {
			print "Less than 4 columns that are at least required, will be skipped: $line\n";
			next;
		}
		
		# Spawn a new process to process the gene
		$pm->start and next;
		
		my $chr = $i[0];
		my $start = $i[1];
		my $end = $i[2];
		my $gene = $i[3];
		if ($start < $extension){
			$start = 0;
		}else{
			$start -= $extension;
		}
		$end += $extension;
		my $plinkname = join(".", $prefix, join("_", $chr, $start, $end, $gene));
		my $plinkfile = File::Spec->catfile($temporary, $plinkname);
		# print "$trait_file $gene $plinkname\n";
		system("plink --bfile $genotype_file --chr $chr --from-bp $start --to-bp $end --maf $maf_filter --indiv-sort 0 --make-bed --out $plinkfile >/dev/null 2>&1");

		if (-e "$plinkfile.bed") {
			# Retry mechanism for failed gene processing
			my ($result_TPED, $result_DNF);
			my $retry_count = 0;
			
			while ($retry_count < $max_retries) {
				($result_TPED, $result_DNF) = function_process_Garfield($plinkfile, $trait_file, $plinkname, $keep_negative);
				$result_TPED =~ s/1\.5 1\.5/1 2/g if $result_TPED=~/1\.5/;;
				last if defined $result_TPED;
				$retry_count++;
			}

			# Store the gene result if successful, otherwise log an error
			if (defined $result_DNF) {
				flock($out_dnf_fh, LOCK_EX);
					print $out_dnf_fh "$result_DNF\n";
				flock($out_dnf_fh, LOCK_UN);
				
				if (defined $out_tped_fh && (not $out_tped_fh=~/NULL/ig)) {
					flock($out_tped_fh, LOCK_EX);
						print $out_tped_fh "$result_TPED\n";
					flock($out_tped_fh, LOCK_UN);
				}
			}
			else {
				warn "Error processing: $ID\n";
			}
		}
		# # # Delete temporary files
		my @plinkfiles = glob("$temporary/$plinkname.*");
		unlink @plinkfiles;
		$pm->finish;
	}
	close $bed_fh;

	$pm->wait_all_children;  # Wait for all child processes to finish
	close $out_tped_fh; close $out_dnf_fh;
}


# Function to process each gene and return the result
sub function_process_Garfield {
	my ($plinkfile, $trait_file, $plinkname, $keep_negative) = @_;
	my $r_output = `Rscript $FindBin::Bin/lib/Garfield_main.functions.R.sh $plinkfile $trait_file $plinkname $keep_negative`;
	my ($result_TPED, $result_DNF) = $r_output =~ /([^\n]+)\n([^\n]+)/;
	return($result_TPED, $result_DNF);
}

1;
