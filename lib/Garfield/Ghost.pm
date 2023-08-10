package Garfield::Ghost;

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
		pod2usage(-exitval => 0, -verbose => 2, -input => "$FindBin::Bin/doc/Garfield_Ghost.txt");
		exit 0;
	}

	if (!check_para($class, $opts_sub)) {
		print "Required parameters missing!\n";
		pod2usage(-exitval => 1, -verbose => 2, -input => "$FindBin::Bin/doc/Garfield_Ghost.txt");
		exit 0;
	}

	my $exit_code = 0;

	if (!$opts_sub->{"INDpeak"}) {
		print "Please provide --INDpeak file!\n";
		$exit_code++;
	} else {
		if (!-e $opts_sub->{"INDpeak"}) {
			print "INDpeak file: $opts_sub->{INDpeak} does not exist. Please check!\n";
			$exit_code++;
		}
		$opts_sub->{"INDpeak"} = abs_path $opts_sub->{"INDpeak"};
	}
	
	
	# Default for Window
	if (!$opts_sub->{"extension"}) {
		$opts_sub->{"extension"} = 100000;
	}
	
	if (!$opts_sub->{"rmLD2peak"}) {
		$opts_sub->{"rmLD2peak"} = 0.3;
	}
	
	if (!$opts_sub->{"LDprune"}) {
		$opts_sub->{"LDprune"} = 0.9;
	}


	#### Common arguments
	my $cm_arg = Garfield::CommonArgument->new();
	my $exit_num_return = $cm_arg->common_argument($opts_sub, $exit_code);
	$exit_code += $exit_num_return;

	if ($exit_code > 0) {
		pod2usage(-exitval => 1, -verbose => 2, -input => "$FindBin::Bin/doc/Garfield_Ghost.txt");
		exit 0;
	} else {
		# return "TRUE";
		print "IND peak file is: $opts_sub->{INDpeak}\n";
		print "With the extension of each INDpeak: $opts_sub->{extension} bp\n";
		print "Variants will be removed if with LD (r2) to the IND peak >=: $opts_sub->{rmLD2peak}\n";
		print "The left variants will be further pruned if with LD (r2) to each other >=: $opts_sub->{LDprune}\n";
		print "\#threads will be used: $opts_sub->{'threads'}\n";
		RUN_Ghost($class, $opts_sub);
	}
}

sub RUN_Ghost {
	my ($class, $opts_sub) = @_;
	my $genotype_file = $opts_sub->{genotype_file};
	# my $trait_file = $opts_sub->{trait_file};
	my $max_threads = $opts_sub->{threads};
	my $INDpeak_file = $opts_sub->{INDpeak};
	my $extension = $opts_sub->{extension};
	my $LD_with_peak = $opts_sub->{rmLD2peak};
	my $LD_prune = $opts_sub->{LDprune};
	my $prefix = $opts_sub->{prefix};
	my $temporary = $opts_sub->{temporary};
	my $output_dir = $opts_sub->{outdir};
	
	open($out_tped_fh, '>', "$output_dir/Garfield.Geno.$prefix.tped") or die "Cannot open $output_dir/Garfield.Geno.$prefix.tped for writing: $!";
	open($out_dnf_fh, '>', "$output_dir/Garfield.bestDNF.$prefix.txt") or die "Cannot open $output_dir/Garfield.bestDNF.$prefix.txt for writing: $!";
	
	my %bim;
	open(my $bim_fh, '<', "$genotype_file.bim") || die "Cannot open $genotype_file.bim: $!";
	while (my $line = <$bim_fh>) {
		chomp $line;
		$line =~ s/\r//g;
		$line =~ s/\ +/\t/g;
		my @i = split("\t", $line);
		if (@i < 6) {
			print "A total of 6 columns are required for plink bim file, please check: $genotype_file.bim\n";
			last;
		}
		my $end = $i[3] + length($i[4]) -1;
		$bim{$i[1]} = join("\t", $i[0], $i[3], $end);
	}
	close $bim_fh;

	# my (%TPED, %bestDNF);
	my $pm = Parallel::ForkManager->new($max_threads);
	
	open(my $INDpeak_fh, '<', $INDpeak_file) || die "Cannot open $INDpeak_file: $!";
	while (my $line = <$INDpeak_fh>) {
		chomp $line;
		$line =~ s/\r//g;
		next if $line =~ /^#/ || $line =~ /^\s*$/;

		if (not exists $bim{$line}) {
			print "This variant not exists in $genotype_file.bim, will be skipped: $line\n";
			next;
		}
		my ($chr, $start, $end) = split("\t", $bim{$line});
		my $snp = $line;
		
		# Spawn a new process to process the gene
		$pm->start and next;

		if ($start < $extension){
			$start = 0;
		}else{
			$start -= $extension;
		}
		$end += $extension;
		
		my $plinkname = join(".", $prefix, join("__", $line));
		my $plinkfile = File::Spec->catfile($temporary, $plinkname);
		# print "$trait_file $gene $plinkname\n";
		system("plink --bfile $genotype_file --chr $chr --from-bp $start --to-bp $end --maf $maf_filter --indiv-sort 0 --make-bed --out $plinkfile >/dev/null 2>&1");
		system("plink --bfile $plinkfile --ld-snp $snp --r2 inter-chr --ld-window-r2 $LD_with_peak --out $plinkfile.LD2peak >/dev/null 2>&1");

		if (-e "$plinkfile.bim" && -e "$plinkfile.LD2peak.ld") {
			open(BIM,"$plinkfile.bim") || die;
			open(LD,"$plinkfile.LD2peak.ld") || die;
			open(OUTLD,">$plinkfile.rmLD.list") || die;
			my %ld;
			$ld{$snp}++;
			while(<LD>){
				next if /CHR_A/;
				s/\ +/\t/g; s/\t+/\t/g; s/^\t//g;
				my @s=split("\t");
				$ld{$s[5]}++;
			}
			while(<BIM>){
				my @s=split("\t");
				print OUTLD "$s[1]\n" unless exists $ld{$s[1]};
			}
			
			system("plink --bfile $plinkfile --extract $plinkfile.rmLD.list --make-bed --out $plinkfile.rmLD >/dev/null 2>&1");
			
			if ($LD_prune<1){
				system("plink --bfile $plinkfile.rmLD --indep-pairwise 10 3 $LD_prune --out $plinkfile >/dev/null 2>&1");
				system("plink --bfile $plinkfile.rmLD --extract $plinkfile.prune.in --make-bed --out $plinkfile.rmLD.prune >/dev/null 2>&1");
			}
			open(OUTPEAK,">$plinkfile.peak") || die;
			print OUTPEAK "$snp\n";
			close OUTPEAK;
			
			system("plink --bfile $plinkfile --extract $plinkfile.peak --recode A --out $plinkfile.peak >/dev/null 2>&1");
			system("cut -d ' ' -f 1,2,7 ${plinkfile}.peak.raw | sed '1d' > ${plinkfile}.trait");

			if ($LD_prune < 1 && -e "$plinkfile.rmLD.prune.bim") {
				process_Ghost_Detail("$plinkfile", 0, $out_tped_fh, $out_dnf_fh);
			}

			if ($LD_prune == 1 && -e "$plinkfile.rmLD.bim") {
				process_Ghost_Detail("$plinkfile", 1, $out_tped_fh, $out_dnf_fh);
			}

		}
		# # Delete temporary files
		my @plinkfiles = glob("$temporary/$plinkname.*");
		unlink @plinkfiles;
		$pm->finish;
	}
	close $INDpeak_fh;

	$pm->wait_all_children;  # Wait for all child processes to finish
	close $out_tped_fh; close $out_dnf_fh;
}


# Function to process each gene and return the result
sub function_process_Garfield {
	my ($plinkfile, $trait_file, $plinkname) = @_;
	my $r_output = `Rscript $FindBin::Bin/lib/Garfield_main.functions.R.sh $plinkfile $trait_file $plinkname`;
	my ($result_TPED, $result_DNF) = $r_output =~ /([^\n]+)\n([^\n]+)/;
	return($result_TPED, $result_DNF);
}


# Common processing for each case
sub process_Ghost_Detail {
    my ($plink_file_prefix, $ld_prune_flag, $out_tped_fh, $out_dnf_fh) = @_;
    
    my $ld_file_suffix = ($ld_prune_flag == 1) ? "rmLD" : "rmLD.prune";
    
    my ($result_TPED, $result_DNF);
    my $retry_count = 0;

    while ($retry_count < $max_retries) {
        ($result_TPED, $result_DNF) = function_process_Garfield("$plink_file_prefix.$ld_file_suffix", "$plinkfile.trait", "$plinkname");
        $result_TPED =~ s/1\.5 1\.5/0 0/g if $result_TPED =~ /1\.5/;
        last if defined $result_TPED;
        $retry_count++;
    }

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



1;
