#!/usr/bin/perl
#
# command to protect CPP options in the rundeck
# when used without options comments out (repalces # with _)
# options between "Preprocessor Options" and "End Preprocessor Options"
#
# when used with -u restores original CPP options
#
# also removes empty lines at the start of rundeck (CPP creates
# extra lines there)

use Getopt::Long;

GetOptions("u") || die "problem in GetOptions";

if ($opt_u ) {
    $unprotect = 1;
}

$start_of_file = 1;

while(<>){

    # skip empy lines at the start of file
    if ($start_of_file && /^\s*$/) { next; }
    $start_of_file = 0;

    if(/^Preprocessor Options/) {
	$inside_cpp_options = 1;
    }
    if(/^End Preprocessor Options/) {
	$inside_cpp_options = 0;
    }

    if( $inside_cpp_options ) {
	if ( (! $unprotect) && /^\#/ ) {
	    s/^\#/_/;
	}
	if ( $unprotect && /^_/ ) {
	    s/^_/\#/;
	}
    }
    print;
}

