#!/usr/bin/perl

if($#ARGV < 0 || $#ARGV > 1) {
    die "usage: $0 C_Source_File [dest]\n";
}

if($#ARGV == 1) {
    $destdir = $ARGV[1];
}
else {
    $destdir = '.';
}

$input = $ARGV[0];
if($input =~ /(.+)\.[^\.]*$/) {
    $basename = $1;
}
else {
    $basename = $input;
}

$header = $destdir . '/' . $basename . '.h';
$headerIN = $basename . '.h.in';

open(header,'>',$header) or die "$header: $!\n";
open(headerIN) or die "$headerIN: $!\n";
while(<headerIN>) {
    last if(/^\s*\@prototypes\s*$/);
    print header;
}

open(input) or die "$input: $!\n";

while(<input>) {
    if(/^(.*)\/\/begin{proto}/) {
        print header "$1\n";
        while(<input>) {
            if(/^(.*)\/\/end{proto}/) {
                $_ = $1;
                s/\s+$//;
                print header "$_;\n\n";
                last;
            }
            print header;
        }
        next;
    }
    if(/^(.*)\/\/prototype/) {
        $_ = $1;
        s/\{//;
        s/\s+$//;
        print header "$_;\n\n";
    }
}
close(input);

while(<headerIN>) {
    last if(/^\s*\@prototypes\s*$/);
    print header;
}

close(headerIN);
close(header);
