package Garfield::Window;

use strict;
use warnings;
use File::Basename;
use FindBin;
use Pod::Usage;
use Cwd qw(abs_path);
use Garfield::CommonArgument;
use Garfield::Process_Window;

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
		return "TRUE";
	}
}

sub run {
    my ($class, $opts_sub) = @_;
	my $cm_arg = Garfield::Process_Window->new();
	$cm_arg->RUN_Window($class, $opts_sub);
}

1;
