#!/usr/bin/perl

use my_bio;
use my_gff;

(@ARGV != 3) and die "usage: $0 genes.gff3 pfam_scan.txt max_e_balue";

$gff3_file = $ARGV[0];
$pfam_file = $ARGV[1];
$max_e_value = $ARGV[2];


gff3_tree_set_relations(\%gff3_tree, $gff3_file);

open($in, $gff3_file);
while (<$in>) {
	chomp;
	($seq_name, $feature, $st, $ed, $strand, $attr) = (split(/\t/, $_))[0, 2, 3, 4, 6, 8];
	next if ($feature ne 'CDS');

	$id = gff3_attr_get_id($attr);
	@mRNA_ids = gff3_tree_get_parent_id(\%gff3_tree, $id);

	for $mRNA_id (@mRNA_ids) {
		$gene{$mRNA_id}{name} = $seq_name;
		$gene{$mRNA_id}{strand} = $strand;
		push(@{$gene{$mRNA_id}{pos}}, $st);
		push(@{$gene{$mRNA_id}{pos}}, $ed);
	}
}
close $in;

open($in, $pfam_file);
while (chomp($l = <$in>)) {
	next if ($l =~ /^#/ or $l eq '');
	($mRNA_id, $type, $e_value) = (split(/\s+/, $l))[0, 7, 12];
	if (($type eq 'Domain' or $type eq 'Family') and $e_value <= $max_e_value) {
		$gene_id = gff3_tree_get_root_id(\%gff3_tree, $mRNA_id);
		$pfam_pass_flag{$gene_id} = 1;
	}
}
close $in;
	

open($in, $gff3_file);
while (<$in>) {
	chomp;
	($attr) = (split(/\t/, $_))[8];

	$id = gff3_attr_get_id($attr);
	$gene_id = gff3_tree_get_root_id(\%gff3_tree, $id);

	if ($pfam_pass_flag{$gene_id}) {
		print "$_\n";
	}
}
close $in;
