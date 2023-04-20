#!/usr/bin/perl

use my_bio;

(@ARGV != 2) and die "usage: $0 gff3 fasta\n";

open IN, $ARGV[0];
while (<IN>) {
	chomp;
	($seq_name, $feature, $st, $ed, $strand, $attr) = (split(/\t/, $_))[0, 2, 3, 4, 6, 8];
	next if ($feature ne 'exon' or $attr =~ /tRNA/);

	$attr =~ /Parent=([^;]+)/;
	$gene_id = $1;
	$gene_id =~ s/\|m\..*$//;
	$gene_id =~ s/-mRNA-\d+$//;

	$gene{$gene_id}{name} = $seq_name;
	$gene{$gene_id}{strand} = $strand;
	push(@{$gene{$gene_id}{pos}}, [$st, $ed]);
}
close IN;

open($in, $ARGV[1]);
while (($name, $s) = fasta_nonwhite_get($in)) {
	$seq{$name} = $s;
}
close $in;


while (($gene_id, $info) = each %gene) {
	@pos = sort {$a->[0] <=> $b->[1]} @{$info->{pos}};
	$mRNA = '';
	for $st_ed (@pos) {
		$mRNA .= substr($seq{$info->{name}}, $st_ed->[0] - 1, $st_ed->[1] - $st_ed->[0] + 1);
	}

	if ($info->{strand} eq '-') {
		$mRNA = rev_comp($mRNA);
	}

	printf(">%s\n", $gene_id);
	print_seq($mRNA);
}
