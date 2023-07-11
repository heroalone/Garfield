package Garfield::Window;

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

my $max_retries = 5;
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
		pod2usage(-exitval => 0, -verbose => 2, -input => "$FindBin::Bin/doc/Garfield_Window.txt");
		exit 0;
	}

	if (!check_para($class, $opts_sub)) {
		print "Required parameters missing!\n";
		pod2usage(-exitval => 1, -verbose => 2, -input => "$FindBin::Bin/doc/Garfield_Window.txt");
		exit 0;
	}

	my $exit_code = 0;

	if (!$opts_sub->{"genome"}) {
		print "Please provide --faidx file!\n";
		$exit_code++;
	} else {
		if (!-e $opts_sub->{"genome"}) {
			print "genome file: $opts_sub->{genome} does not exist. Please check!\n";
			$exit_code++;
		}
		$opts_sub->{"genome"} = abs_path $opts_sub->{"genome"};
	}

	# Default for Window
	if (!$opts_sub->{"window"}) {
		$opts_sub->{"window"} = 50000;
	}
	
	if (!$opts_sub->{"step"}) {
		$opts_sub->{"step"} = 25000;
	}

	#### Common arguments
	my $cm_arg = Garfield::CommonArgument->new();
	my $exit_num_return = $cm_arg->common_argument($opts_sub, $exit_code);
	$exit_code += $exit_num_return;

	if ($exit_code > 0) {
		pod2usage(-exitval => 1, -verbose => 2, -input => "$FindBin::Bin/doc/Garfield_Window.txt");
		exit 0;
	} else {
		# return "TRUE";
		RUN_Window($class, $opts_sub);
	}
}


sub RUN_Window{
	my ($class, $opts_sub) = @_;
	my $genotype_file = $opts_sub->{genotype_file};
	my $trait_file = $opts_sub->{trait_file};
	my $max_threads = $opts_sub->{threads};
	my $genome_file = $opts_sub->{genome};
	my $window = $opts_sub->{window};
	my $step = $opts_sub->{step};
	my $temporary = $opts_sub->{temporary};
	# my (%TPED, %bestDNF);
	my $pm;
	my $output_dir = $opts_sub->{outdir};
	my $prefix = $opts_sub->{prefix};
	
	open($out_tped_fh, '>', "$output_dir/Garfield.Geno.$prefix.tped") or die "Cannot open $output_dir/Garfield.Geno.$prefix.tped for writing: $!";
	open($out_dnf_fh, '>', "$output_dir/Garfield.bestDNF.$prefix.txt") or die "Cannot open $output_dir/Garfield.bestDNF.$prefix.txt for writing: $!";
	
	open(my $genome_fh, '<', $genome_file) || die "Cannot open $genome_file: $!";
	while (my $line = <$genome_fh>) {
		chomp $line;
		$line =~ s/\r//g;
		next if $line =~ /^#/ || $line =~ /^\s*$/;
		
		my ($chr, $genome_length, $tmp) = split("\t", $line, 3);
		
		for (my $start = 0; $start < $genome_length; $start += $step) {
			my $end = $start + $window;
			# Spawn a new process to process the Window
			
			$pm = Parallel::ForkManager->new($max_threads);
			$pm->start and next;

			my $ID = join("_", $chr, $start, $end);
			my $plinkname=join(".",$prefix, $ID);
			my $plinkfile = File::Spec->catfile($temporary, $plinkname);
			# # print "$trait_file $plinkfile\n";

			system("plink --bfile $genotype_file --chr $chr --from-bp $start --to-bp $end --maf $maf_filter --make-bed --out $plinkfile >/dev/null 2>&1");
			if (-e "$plinkfile.bed") {
				# Retry mechanism for failed Window processing
				my ($result_TPED, $result_DNF);
				my $retry_count = 0;
				
				while ($retry_count < $max_retries) {
					($result_TPED, $result_DNF) = function_process_Garfield($plinkfile, $trait_file, $plinkname);
					last if defined $result_TPED;
					$retry_count++;
				}
				
				# Store the Window result if successful, otherwise log an error
				if (defined $result_TPED) {
					flock($out_tped_fh, LOCK_EX);
					flock($out_dnf_fh, LOCK_EX);
						print $out_tped_fh "$result_TPED\n";
						print $out_dnf_fh "$result_DNF\n";
					flock($out_tped_fh, LOCK_UN);
					flock($out_dnf_fh, LOCK_UN);
				} else {
					warn "Error processing: $ID\n";
				}
			}
			
			# # # Delete temporary files
			my @plinkfiles = glob("$temporary/$plinkname.*");
			unlink @plinkfiles;
			
			$pm->finish;
		}
	}
	close $genome_fh;
	
	$pm->wait_all_children;  # Wait for all child processes to finish
	close $out_tped_fh; close $out_dnf_fh;
}


# Function to process each Window and return the result
sub function_process_Garfield {
	my ($plinkfile, $trait_file, $plinkname) = @_;
	my $r_output = `Rscript $FindBin::Bin/lib/Garfield_main.functions.R.sh $plinkfile $trait_file $plinkname`;
	my ($result_TPED, $result_DNF) = $r_output =~ /([^\n]+)\n([^\n]+)/;
	return($result_TPED, $result_DNF);
}

1;
