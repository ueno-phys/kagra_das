#!/usr/bin/perl

if($#ARGV != 0) {
    die "usage: $0 name\n";
}

$basename = $ARGV[0];
$output = $basename . 'Version.h';

open(OUT,'>',$output) or die "output: $!\n";

$gitver = `/usr/bin/git log -1 --format="%ci %cn %H"`;
chomp $gitver;
$gitid = `/usr/bin/git log -1 --format="%h"`;
chomp $gitid;
print OUT <<__END__;
#define KGL_GITVER "\$$basename: $gitver \$"
#define KGL_GITID "$gitid"
__END__
close(OUT);
