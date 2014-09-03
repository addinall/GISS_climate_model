#!/usr/bin/perl

# to skip dependencies on certan modules list them here
# for example:
# $skip_modules="domain_decomp.mod|model_com.mod|constant.mod";
$skip_modules="";

while(<>) {

    # hack to skip unneeded dependencies
    if ( $skip_modules ) { s/ ($skip_modules)//g; }

    $MOD = "mod"; # could be different for other compilers

    if ( /\s*(\w+)\.$MOD:\s+(\w+)\.o\s*$/ ) {
	print "$1\@$2.smod: $2.o\n";
	print "$1.$MOD: $1\@$2.smod\n";
    } else {
	print ;
    }


}
