package Garfield::CommonArgument;

use strict;
use warnings;
use File::Spec;
use File::Path;
use Cwd qw(abs_path);

sub new {
	my $class = shift;
	return bless {}, $class;
}

sub common_argument {
	my ($class, $opts_sub, $exit_code) = @_;
	my ($sub_cmd) = @{$opts_sub->{subcommand}};
	my ($tmp_geno, $tmp_pheno);

	if (!$opts_sub->{"genotype_file"}) {
		print "Please provide --genotype file!\n";
		$exit_code++;
	} else {
		$tmp_geno=$opts_sub->{"genotype_file"};
		
		if (!-e "$tmp_geno\.bed") {
			print "genotype file: $tmp_geno.bed does not exist. Please check!\n";
			$exit_code++;
		}elsif (!-e "$tmp_geno\.bim") {
			print "genotype file: $tmp_geno.bim does not exist. Please check!\n";
			$exit_code++;
		}elsif (!-e "$tmp_geno\.fam") {
			print "genotype file: $tmp_geno.fam does not exist. Please check!\n";
			$exit_code++;
		}
		$opts_sub->{"genotype_file"} = abs_path $opts_sub->{"genotype_file"};
	}

	if (!$opts_sub->{"trait_file"}) {
		print "Please provide --trait file!\n";
		$exit_code++;
	} else {
		$tmp_pheno=$opts_sub->{"trait_file"};
		if (!-e $opts_sub->{"trait_file"}) {
			print "phenotype file: $opts_sub->{trait_file} does not exist. Please check!\n";
			$exit_code++;
		}
		$opts_sub->{"trait_file"} = abs_path $opts_sub->{"trait_file"};
	}
	
	#####################
	if (!$opts_sub->{"outdir"}) {
		$opts_sub->{"outdir"} = abs_path "./";
	} else {
		if (-e $opts_sub->{"outdir"} && !-d $opts_sub->{"outdir"}) {
			print "Directory $opts_sub->{outdir} already exists. Please provide a new directory name.\n";
			$exit_code++;
		}
		unless (-d $opts_sub->{"outdir"}) {
			File::Path::make_path($opts_sub->{"outdir"}) or die "Failed to create directory: $!";
		}
		$opts_sub->{"outdir"} = abs_path $opts_sub->{"outdir"};
	}
	
	unless (-d $opts_sub->{"outdir"}) {
		File::Path::make_path($opts_sub->{"outdir"}) or die "Failed to create directory: $!";
	}
	print "Output directory is: $opts_sub->{outdir}\n" if $exit_code == 0;


	#####################
	if (!$opts_sub->{"temporary"}) {
		$opts_sub->{"temporary"} = abs_path "./tmp";
	} else {
		if (-e $opts_sub->{"temporary"} && !-d $opts_sub->{"temporary"}) {
			print "Temporary directory $opts_sub->{temporary} already exists. Please provide a new directory name.\n";
			$exit_code++;
		}
		unless (-d $opts_sub->{"temporary"}) {
			File::Path::make_path($opts_sub->{"temporary"}) or die "Failed to create temporary directory: $!";
		}
		$opts_sub->{"temporary"} = abs_path $opts_sub->{"temporary"};
	}
	
	unless (-d $opts_sub->{"temporary"}) {
		File::Path::make_path($opts_sub->{"temporary"}) or die "Failed to create temporary directory: $!";
	}
	print "Temporary directory is: $opts_sub->{temporary}\n" if $exit_code == 0;
	#####################


	if ($opts_sub->{"genotype_file"} && $opts_sub->{"trait_file"}) {
		if (!$opts_sub->{"prefix"}) {
			$opts_sub->{"prefix"} = join("_",  $tmp_geno,  $tmp_pheno, $sub_cmd);
			print "Default output prefix is used: $opts_sub->{'prefix'}\n" if $exit_code == 0;
		} else {
			print "Output prefix: $opts_sub->{'prefix'}\n" if $exit_code == 0;
		}
	}

	if (!$opts_sub->{"threads"}) {
		$opts_sub->{"threads"} = 1;
	}
	return $exit_code;
}

1;
