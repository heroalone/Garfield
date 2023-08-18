#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use Cwd qw(abs_path);

# List of Perl modules to be installed
my @perl_modules = ("Benchmark", "FindBin", "File::Spec", "Pod::Usage", "Getopt::Long::Subcommand", "Fcntl", "Parallel::ForkManager");
my $dir = dirname(abs_path $0);


# Get the path of cpanm
my $cpanm_path;
$cpanm_path = qx(which cpanm 2>/dev/null);

if ($? == 0 && $cpanm_path) {
	print "use cpanm at: $cpanm_path\n";
}
else {
	my $script_cpanm = "$dir/scripts/cpanminus";
	if (!-e $script_cpanm) {
		die "Error: cpanm not found or not executable in your system path and the accompanying package. Please install cpanm first to proceed.\n";
	}
	else {
		$cpanm_path = $script_cpanm;
		system("chmod +x $cpanm_path");
		print "use cpanm at: $cpanm_path\n";
	}
}


# Check if local::lib is set up
my $perl_lib_check = qx(perl -e 'use local::lib;');
if (!$perl_lib_check) {
	system("$cpanm_path --local-lib=~/perl5 local::lib && eval \$(perl -I ~/perl5/lib/perl5/ -Mlocal::lib)");
}

# Install Perl modules
foreach my $module (@perl_modules) {
	eval "use $module;";
	if ($@) {
		print "Installing Perl module $module...\n";
		system("$cpanm_path $module");
		if ($? == 0) {
			print "Perl module $module installed successfully.\n";
		} else {
			die "Failed to install Perl module $module. Please install it manually.\n";
		}
	}
	else {
		print "Perl module $module is already installed.\n";
	}
}

# Check for R packages
my $rscript_path = checkINpath("Rscript");
if (!$rscript_path) {
	die "Rscript not found in the PATH. Please install R or add it to your PATH.\n";
}
else {
	print "Rscript found in the PATH.\n";
}

# Run R package installation
my $r_package_script = "$dir/scripts/install_packages.R.sh";
my $r_package_output = qx(Rscript $r_package_script 2>&1);
print "$r_package_output\n";

if ($r_package_output =~ /successful/) {
    print "All dependencies have been satisfied.\n";
    print "Now you can take Garfield on your journey!\n";
	
	# my $ascii = "$dir/images/Garfield.ascii_art";
    if (-e "$dir/images/Garfield.ascii_art") {
		open my $ascii_file, '<', "$dir/images/Garfield.ascii_art";
		while (my $line = <$ascii_file>) {
			print $line;
		}
		close $ascii_file;
	}
}



sub checkINpath {
	my ($tool_name) = @_;
	my @path_dirs = split /:/, $ENV{PATH};
	foreach my $path (@path_dirs) {
		if (-x "$path/$tool_name") {
			return $path;
		}
	}
	return undef;
}
