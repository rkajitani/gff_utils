#!/usr/bin/perl

use my_bio;

(@ARGV != 2) and die "usage: $0 gff3 fasta\n";

open IN, $ARGV[0];
while (<IN>) {
	chomp;
	($seq_name, $feature, $st, $ed, $strand, $attr) = (split(/\t/, $_))[0, 2, 3, 4, 6, 8];
	next if ($feature ne 'CDS' or $attr =~ /tRNA/);

	$attr =~ /ID=([^;]*)/;
	$gene_id = $1;

	$gene{$gene_id}{name} = $seq_name;
	$gene{$gene_id}{strand} = $strand;
	push(@{$gene{$gene_id}{pos}}, $st);
	push(@{$gene{$gene_id}{pos}}, $ed);
}
close IN;

open($in, $ARGV[1]);
while (($name, $s) = fasta_nonwhite_get($in)) {
	$seq{$name} = $s;
}
close $in;


while (($gene_id, $info) = each %gene) {
	@pos = sort {$a <=> $b} @{$info->{pos}};
	$CDS = '';
	for ($i = 0; $i < $#pos; $i += 2) {
		$CDS .= substr($seq{$info->{name}}, $pos[$i] - 1, $pos[$i + 1] - $pos[$i] + 1);
	}

	if ($info->{strand} eq '-') {
		$CDS = rev_comp($CDS);
	}

	printf(">%s\n", $gene_id);
	print_seq($CDS);
}
