package Garfield::GeneSet;

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
		pod2usage(-exitval => 0, -verbose => 2, -input => "$FindBin::Bin/doc/Garfield_GeneSet.txt");
		exit 0;
	}

	if (!check_para($class, $opts_sub)) {
		print "Required parameters missing!\n";
		pod2usage(-exitval => 1, -verbose => 2, -input => "$FindBin::Bin/doc/Garfield_GeneSet.txt");
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
	
	if (!$opts_sub->{"geneset"}) {
		print "Please provide --geneset file!\n";
		$exit_code++;
	} else {
		if (!-e $opts_sub->{"geneset"}) {
			print "geneset file: $opts_sub->{geneset} does not exist. Please check!\n";
			$exit_code++;
		}
		$opts_sub->{"geneset"} = abs_path $opts_sub->{"geneset"};
	}
	

	# Default for GeneSet
	if (!$opts_sub->{"extension"}) {
		$opts_sub->{"extension"} = 20000;
	}

	#### Common arguments
	my $cm_arg = Garfield::CommonArgument->new();
	my $exit_num_return = $cm_arg->common_argument($opts_sub, $exit_code);
	$exit_code += $exit_num_return;

	if ($exit_code > 0) {
		pod2usage(-exitval => 1, -verbose => 2, -input => "$FindBin::Bin/doc/Garfield_GeneSet.txt");
		exit 0;
	} else {
		# return "TRUE";
		RUN_GeneSet($class, $opts_sub);
	}
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
	my $output_dir = $opts_sub->{outdir};
	
	open($out_tped_fh, '>', "$output_dir/Garfield.Geno.$prefix.tped") or die "Cannot open $output_dir/Garfield.Geno.$prefix.tped for writing: $!";
	open($out_dnf_fh, '>', "$output_dir/Garfield.bestDNF.$prefix.txt") or die "Cannot open $output_dir/Garfield.bestDNF.$prefix.txt for writing: $!";
	
	
	my %bed; my $pm;
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
		$pm = Parallel::ForkManager->new($max_threads);
		$pm->start and next;


		my $plinkname = join(".", $prefix, join("__", @i));
		my $plinkfile = File::Spec->catfile($temporary, $plinkname);
		for (my $n=0; $n<@i; $n++){
				unless (exists $bed{$i[$n]}){
					print "ERROR: $i[$n] is present in GeneSet but not in gene bed file, please check!\n";
					exit;
				}
				my ($chr, $start, $end)=split("\t", $bed{$i[$n]});
				
				$start -= $extension if $start >= $extension;
				$start=0 if $start < $extension;
				$end += $extension;
				system("plink --bfile $genotype_file --chr $chr --from-bp $start --to-bp $end --maf $maf_filter --write-snplist --out $plinkfile.$i[$n] >/dev/null 2>&1");
				system("cat $plinkfile.$i[$n].snplist >> $temporary/merge.$plinkname.snplist");
				system ("rm $plinkfile.$i[$n].*");
			}
			system("plink --bfile $genotype_file --extract $temporary/merge.$plinkname.snplist --make-bed --out $plinkfile >/dev/null 2>&1");
		if (-e "$plinkfile.bed") {
			my ($result_TPED, $result_DNF);
			my $retry_count = 0;
			while ($retry_count < $max_retries) {
				my $r_output = `Rscript $FindBin::Bin/lib/Garfield_main.functions.R.sh $plinkfile $trait_file $plinkname`;
				
				($result_TPED, $result_DNF) =  $r_output =~ /([^\n]+)\n([^\n]+)/;
				last if defined $result_TPED;
				$retry_count++;
			}

			# Store the geneset result if successful, otherwise log an error
			if (defined $result_TPED) {
				flock($out_tped_fh, LOCK_EX);
				flock($out_dnf_fh, LOCK_EX);
					print $out_tped_fh "$result_TPED\n";
					print $out_dnf_fh "$result_DNF\n";
				flock($out_tped_fh, LOCK_UN);
				flock($out_dnf_fh, LOCK_UN);
			} else {
				warn "Error processing: $line\n";
			}
		}
		# # # Delete temporary files
		my @plinkfiles = glob("$temporary/$plinkname.*");
		unlink @plinkfiles;
		
		$pm->finish;

	}
	close $geneset_fh;

	$pm->wait_all_children; 
	close $out_tped_fh;close $out_dnf_fh;
}

1;
