package Garfield::Process_GeneSet;

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

sub RUN_GeneSet{
	my ($class, $opts_sub) = @_;
	my $genotype_file = $opts_sub->{genotype_file};
	my $trait_file = $opts_sub->{trait_file};
	my $max_threads = $opts_sub->{threads};
	my $bed_file = $opts_sub->{bed};
	my $extension = $opts_sub->{extension};
	my $prefix = $opts_sub->{prefix};
	my $geneset_file = $opts_sub->{geneset};
	my $temporary = $opts_sub->{temporary};
	
	my (%TPED, %bestDNF, %bed);
	

	open(my $bed_fh, '<', $bed_file) || die "Cannot open $bed_file: $!";
	while (my $line = <$bed_fh>) {
		chomp $line;
		$line =~ s/\r//g;
		next if $line =~ /^#/ || $line =~ /^\s*$/;

		my @i = split("\t", $line);
		if (@i < 4) {
			print "Less than 4 columns that are at least required, will be missed: $line\n";
			next;
		}
		$bed{$i[3]} = join("\t", $i[0], $i[1], $i[2]);
	}
	close $bed_fh;
	
	open(my $geneset_fh, '<', $geneset_file) || die "Cannot open $geneset_file: $!";
	while (my $line = <$geneset_fh>) {
		chomp $line;
		$line =~ s/\r//g;
		next if $line =~ /^#/ || $line =~ /^\s*$/;

		my @i = split("\t", $line);
		if (@i < 2) {
			print "At least 2 genes should be input, will be missed: $line\n";
			next;
		}
		
		# Spawn a new process to process the geneset
		my $pm = Parallel::ForkManager->new($max_threads);
		$pm->start and next;

		my ($result_TPED, $result_DNF);
		my $retry_count = 0;
		my $plinkname;
		while ($retry_count < $max_retries) {
			$plinkname=join(".", $prefix, join("__", @i));
			my $plinkfile = File::Spec->catfile($temporary, $plinkname);
			# print "$TRAIT $BED $plinkname\n";

			for (my $n=0; $n<@i; $n++){
				unless (exists $bed{$i[$n]}){
					print "ERROR: $i[$n] is present in GeneSet but not in gene bed file, please check!\n";
					exit;
				}
				my ($chr, $start, $end)=split("\t", $bed{$i[$n]});
				
				$start -= $extension if $start >= $extension;
				$start=0 if $start < $extension;
				$end += $extension;
				system("plink --bfile $genotype_file --chr $chr --from-bp $start --to-bp $end --maf $maf_filter --write-snplist --out $plinkfile.$i[$n]");
				system("cat $plinkfile.$i[$n].snplist >> $temporary/merge.$plinkname.snplist");
				system ("rm $plinkfile.$i[$n].*");
			}
			system("plink --bfile $genotype_file --extract $temporary/merge.$plinkname.snplist --make-bed --out $plinkfile");
			my $r_output = `Rscript $FindBin::Bin/lib/Garfield_main.functions.R.sh $plinkfile $trait_file $plinkname`;
			
			($result_TPED, $result_DNF) =  $r_output =~ /([^\n]+)\n([^\n]+)/;
			
			my @plinkfiles = glob("$temporary/$plinkfile.*");
			unlink @plinkfiles;
			
			last if defined $result_TPED;
			$retry_count++;
		}

		# Store the geneset result if successful, otherwise log an error
		if (defined $result_TPED) {
			$TPED{$plinkname} = $result_TPED;
			$bestDNF{$plinkname} = $result_DNF;

		} else {
			warn "Error processing: $line\n";
		}

		$pm->finish;

	}
	close $geneset_fh;

	my $output_dir = $opts_sub->{outdir};
	my $output_prefix = $opts_sub->{prefix};
	open(my $out_tped_fh, '>', "$output_dir/Garfield.Geno.$output_prefix.tped") or die "Cannot open $output_dir/Garfield.Geno.$output_prefix.tped for writing: $!";
	open(my $out_dnf_fh, '>', "$output_dir/Garfield.bestDNF.$output_prefix.txt") or die "Cannot open $output_dir/Garfield.bestDNF.$output_prefix.txt for writing: $!";

	# Write the geneset results to the output file
	foreach my $plinkname9 (sort keys %TPED) {
		print $out_tped_fh "$TPED{$plinkname9}\n";
		print $out_dnf_fh "$bestDNF{$plinkname9}\n";
	}
	close $out_tped_fh;close $out_dnf_fh;

}

1;