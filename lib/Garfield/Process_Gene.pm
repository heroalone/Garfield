package Garfield::Process_Gene;

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

sub RUN_Gene{
	my ($class, $opts_sub) = @_;
	# my $genotype_file = $opts_sub->{genotype_file};
	# my $trait_file = $opts_sub->{trait_file};
	my $max_threads = $opts_sub->{threads};
	my $bed_file = $opts_sub->{bed};
	# my $extension = $opts_sub->{extension};
	# my $prefix = $opts_sub->{prefix};
	my (%TPED, %bestDNF, $plinkname);
	

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
		my $pm = Parallel::ForkManager->new($max_threads);
		$pm->start and next;
		
		$plinkname = join(".", $opts_sub->{prefix}, join("_", $i[0], $i[1], $i[2], $i[3]));
		
		# Retry mechanism for failed gene processing
		my ($result_TPED, $result_DNF);
		my $retry_count = 0;
		
		while ($retry_count < $max_retries) {
			($result_TPED, $result_DNF) = function_process_gene($i[0], $i[1], $i[2], $i[3], $opts_sub->{extension}, $opts_sub->{prefix}, $opts_sub->{genotype_file}, $opts_sub->{trait_file}, $opts_sub->{temporary});
			
			last if defined $result_TPED;
			$retry_count++;
		}

		# Store the gene result if successful, otherwise log an error
		if (defined $result_TPED) {
			$TPED{$plinkname} = $result_TPED;
			$bestDNF{$plinkname} = $result_DNF;

		} else {
			warn "Error processing: $_\n";
		}
		$pm->finish;

	}
	close $bed_fh;

	my $output_dir = $opts_sub->{outdir};
	my $output_prefix = $opts_sub->{prefix};
	my ($out_tped_fh, $out_dnf_fh);
	open($out_tped_fh, '>', "$output_dir/Garfield.Geno.$output_prefix.tped") or die "Cannot open $output_dir/Garfield.Geno.$output_prefix.tped for writing: $!";
	open($out_dnf_fh, '>', "$output_dir/Garfield.bestDNF.$output_prefix.txt") or die "Cannot open $output_dir/Garfield.bestDNF.$output_prefix.txt for writing: $!";

	# Write the gene results to the output file
	foreach my $plinkname (sort keys %TPED) {
		print $out_tped_fh "$TPED{$plinkname}\n";
		print $out_dnf_fh "$bestDNF{$plinkname}\n";
	}
	close $out_tped_fh;
	close $out_dnf_fh;

}


# Function to process each gene and return the result
sub function_process_gene {
	my ($chr, $start, $end, $gene, $extension, $prefix, $genotype_file, $trait_file, $temporary) = @_;
	if ($start < $extension){
		$start = 0;
	}else{
		$start -= $extension;
	}
	$end += $extension;
	my $plinkname = join(".", $prefix, join("_", $chr, $start, $end, $gene));
	my $plinkfile = File::Spec->catfile($temporary, $plinkname);
	print "$trait_file $gene $plinkname\n";


	system("plink --bfile $genotype_file --chr $chr --from-bp $start --to-bp $end --maf $maf_filter --make-bed --out $plinkfile");
	my $r_output = `Rscript $FindBin::Bin/lib/Garfield_main.functions.R.sh $plinkfile $trait_file $plinkname`;
	my ($result_TPED, $result_DNF) = $r_output =~ /([^\n]+)\n([^\n]+)/;
	

	# Delete temporary files
	my @plinkfiles = glob("$temporary/$plinkname.*");
	unlink @plinkfiles;
	
	return($result_TPED, $result_DNF);
}

1;