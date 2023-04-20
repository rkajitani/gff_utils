#!/usr/bin/perl

use my_bio;
use my_gff;

(@ARGV != 5) and die "usage: $0 genes.gff3 blast.outfmt6 query.fa max_e_value min_query_length";

$max_e_value = $ARGV[3];
$min_query_len = $ARGV[4];

gff3_tree_set_relations(\%gff3_tree, $ARGV[0]);

open($in, $ARGV[1]);
while (chomp($l = <$in>)) {
	($q_name, $t_name, $e_value) = (split(/\t/, $l))[0, 1, 10];	

	if ($e_value <= $max_e_value) {
		$hit_mRNA_flag{$q_name} = 1;
	}
}
close $in;

open($in, $ARGV[2]);
while (($name, $seq) = fasta_nonwhite_get($in)) {
	$seq_len{$name} = length($seq);
}
close $in;

open($in, $ARGV[0]);
while (<$in>) {
	chomp;
	($seq_name, $feature, $st, $ed, $strand, $attr) = (split(/\t/, $_))[0, 2, 3, 4, 6, 8];
	next if ($feature ne 'mRNA');

	$id = gff3_attr_get_id($attr);
	if (not($hit_mRNA_flag{$id}) and $seq_len{$id} < $min_query_len) {
		$gene_id = gff3_tree_get_root_id(\%gff3_tree, $id);
		$removal_flag{$gene_id} = 1;
	}
}
close $in;

open($in, $ARGV[1]);
while (($name, $s) = fasta_nonwhite_get($in)) {
	$seq{$name} = $s;
}
close $in;

open($in, $ARGV[0]);
while (<$in>) {
	chomp;
	($attr) = (split(/\t/, $_))[8];

	$id = gff3_attr_get_id($attr);
	$gene_id = gff3_tree_get_root_id(\%gff3_tree, $id);

	if (not $removal_flag{$gene_id}) {
		print "$_\n";
	}
}
close $in;
