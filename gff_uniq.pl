#!/usr/bin/perl

use my_bio;

(@ARGV != 1) and die "usage: $0 gene.gff3\n";

open IN, $ARGV[0];
while (chomp($l = <IN>)) {
	next if ($l =~ /#/);
	($feature, $st, $ed, $attr) = (split(/\t/, $l))[2, 3, 4, 8];
	if ($feature eq 'mRNA') {
		$attr =~ /ID=([^;]*)/;
		$id = $1;
		$attr =~ /Parent=([^;]*)/;
		$parent = $1;
		if ($ed - $st > $longest{$parent}) {
			 $longest{$parent} = $ed - $st;
			 $is_longest_mRNA{$longest_mRNA{$parent}} = 0;;
			 $is_longest_mRNA{$id} = 1;
			 $longest_mRNA{$parent} = $id;
		}
	}
}
close IN;

open IN, $ARGV[0];
while (chomp($l = <IN>)) {
	($feature, $st, $ed, $attr) = (split(/\t/, $l))[2, 3, 4, 8];
	$attr =~ /ID=([^;]*)/;
	$id = $1;
	$attr =~ /Parent=([^;]*)/;
	$parent = $1;

	if ($feature eq 'gene' or $l =~ /#/) {
		print "$l\n";
	}
	elsif ($feature eq 'mRNA') {
		print "$l\n", if ($is_longest_mRNA{$id});
	}
	elsif ($parent =~ /,/) {
		@parents = split(/,/, $parent);
		for (@parents) {
			if ($is_longest_mRNA{$_}) {
				print "$l\n";
				last;
			}
		}
	}
	elsif (not $is_longest_mRNA{$parent}) {
		next;
	}
	else {
		print "$l\n";
	}
}
close IN;
