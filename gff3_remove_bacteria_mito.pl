#!/usr/bin/perl

(@ARGV != 1) and die "usage: $0 gff3\n";

open($in, $ARGV[0]);
while (chomp($l = <$in>)) {
	($seq_name) = split(/\t/, $l);
	if ($seq_name !~ /Bacteria/ and $seq_name !~ /Mitochondrion/) {
		print "$l\n";
	}
}
close $in;
