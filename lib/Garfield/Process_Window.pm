package Garfield::Process_Window;

use strict;
use warnings;
use Parallel::ForkManager;
use File::Spec;

my $max_retries = 5;
my $maf_filter = 0.05;

sub new {
	my $class = shift;
	return bless {}, $class;
}

sub RUN_Window{
	my ($class, $opts_sub) = @_;
	# my $genotype_file = $opts_sub->{genotype_file};
	# my $trait_file = $opts_sub->{trait_file};
	my $max_threads = $opts_sub->{threads};
	my $genome_file = $opts_sub->{genome};
	my $window = $opts_sub->{window};
	my $step = $opts_sub->{step};
	my (%TPED, %bestDNF);
	

	open(my $genome_fh, '<', $genome_file) || die "Cannot open $genome_file: $!";
	while (my $line = <$genome_fh>) {
		chomp $line;
		$line =~ s/\r//g;
		next if $line =~ /^#/ || $line =~ /^\s*$/;
		
		my ($chr, $genome_length, $tmp) = split("\t", $line, 3);
		my $start = 0;
		
		while ($start + $window <= $genome_length) {
			my $end = $start + $window;
			my $id = join("_", $chr, $start, $end);
			my $plinkname = join(".", $opts_sub->{prefix}, $id);
		
			# Spawn a new process to process the Window
			my $pm = Parallel::ForkManager->new($max_threads);
			$pm->start and next;

			# Retry mechanism for failed Window processing
			my ($result_TPED, $result_DNF);
			my $retry_count = 0;
			
			while ($retry_count < $max_retries) {
				($result_TPED, $result_DNF) = function_process_Window($chr, $start, $end, $id, $opts_sub->{prefix}, $opts_sub->{genotype_file}, $opts_sub->{trait_file}, $opts_sub->{temporary});
				last if defined $result_TPED;
				$retry_count++;
			}

			# Store the Window result if successful, otherwise log an error
			if (defined $result_TPED) {
				$TPED{$plinkname} = $result_TPED;
				$bestDNF{$plinkname} = $result_DNF;
			} else {
				warn "Error processing: $id\n";
			}

			
			$start += $step;
			$pm->finish;
		}

	}
	close $genome_fh;

	my $output_dir = $opts_sub->{outdir};
	my $output_prefix = $opts_sub->{prefix};
	open(my $out_tped_fh, '>', "$output_dir/Garfield.Geno.$output_prefix.tped") or die "Cannot open $output_dir/Garfield.Geno.$output_prefix.tped for writing: $!";
	open(my $out_dnf_fh, '>', "$output_dir/Garfield.bestDNF.$output_prefix.txt") or die "Cannot open $output_dir/Garfield.bestDNF.$output_prefix.txt for writing: $!";

	# Write the Window results to the output file
	foreach my $plinkname (sort keys %TPED) {
		print $out_tped_fh "$TPED{$plinkname}\n";
		print $out_dnf_fh "$bestDNF{$plinkname}\n";
	}
	close $out_tped_fh; close $out_dnf_fh;

}


# Function to process each Window and return the result
sub function_process_Window {
	my ($chr, $start, $end, $ID, $prefix, $genotype_file, $trait_file, $temporary) = @_;

	my $plinkname=join(".",$prefix, $ID);
	my $plinkfile = File::Spec->catfile($temporary, $plinkname);
	print "$trait_file $plinkfile\n";
	

	system("plink --bfile $genotype_file --chr $chr --from-bp $start --to-bp $end --maf $maf_filter --make-bed --out $plinkfile");
	
	my $r_output = `Rscript $FindBin::Bin/lib/Garfield_main.functions.R.sh $plinkfile $trait_file $plinkname`;

	my ($result_TPED, $result_DNF) = $r_output =~ /([^\n]+)\n([^\n]+)/;
	
	my @plinkfiles = glob("$temporary/$plinkname.*");
	unlink @plinkfiles;
	
	return($result_TPED, $result_DNF);
}

1;