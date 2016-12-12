#!/usr/bin/env perl
use strict;
use warnings;
use constant USAGE =><<EOH;

usage: $0 input.synteny 1000 output.synteny

v20161129

EOH
die USAGE if (scalar(@ARGV) !=3 or $ARGV[0] eq '-h' or $ARGV[0] eq '--help');



open (SYNTENY, "< $ARGV[0]") || die "Error: can not open input\n";
my $maxdist=$ARGV[1];
open (OUTPUT, " > $ARGV[2]") || die "Error: can not write output\n";

my $lastqueryid='';
my $lastquerystart='';
my $lastqueryend='';
my $lastsubjectid='';
my $lastsubjectstart='';
my $lastsubjectend='';
my $laststrand='';
my $linenum=0;

my @syntenyarr=();
while (my $line=<SYNTENY>) {
	$linenum++;
	chomp $line;
	next if ($line=~/^#/);
	@syntenyarr=();
	@syntenyarr=split(/\t/, $line);
	if ($lastqueryid eq '') {
		$lastqueryid=$syntenyarr[0]; $lastquerystart=$syntenyarr[1]; $lastqueryend=$syntenyarr[2];
		$lastsubjectid=$syntenyarr[3]; $lastsubjectstart=$syntenyarr[4]; $lastsubjectend=$syntenyarr[5];
		$laststrand=$syntenyarr[7];
		next;
	}
	if ($lastqueryid eq $syntenyarr[0] and $lastsubjectid eq $syntenyarr[3] and $laststrand eq $syntenyarr[7]) {
		if (($syntenyarr[1]>=$lastqueryend) and (($syntenyarr[1]-$lastqueryend)<=$maxdist) and ($syntenyarr[4]>=$lastsubjectend) and (($syntenyarr[4]-$lastsubjectend)<=$maxdist)) {
			$lastqueryend=$syntenyarr[2];
			$lastsubjectend=$syntenyarr[5];
			next;
		}
		else {
			&PrintLast;
		}
	}
	else {
		&PrintLast;
	}

}
print OUTPUT $lastqueryid, "\t", $lastquerystart, "\t", $lastqueryend, "\t", $lastsubjectid, "\t", $lastsubjectstart, "\t", $lastsubjectend, "\t", ($lastqueryend-$lastquerystart+1), "\t", $laststrand;
close SYNTENY;
close OUTPUT;


sub PrintLast {

	print OUTPUT $lastqueryid, "\t", $lastquerystart, "\t", $lastqueryend, "\t", $lastsubjectid, "\t", $lastsubjectstart, "\t", $lastsubjectend, "\t", ($lastqueryend-$lastquerystart+1), "\t", $laststrand, "\n";
	$lastqueryid=$syntenyarr[0]; $lastquerystart=$syntenyarr[1]; $lastqueryend=$syntenyarr[2];
	$lastsubjectid=$syntenyarr[3]; $lastsubjectstart=$syntenyarr[4]; $lastsubjectend=$syntenyarr[5];
	$laststrand=$syntenyarr[7];
}
