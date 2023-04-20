#!/usr/bin/perl

use my_bio;
use my_gff;

(@ARGV != 5) and die "usage: $0 genes.gff3 genome.fa pfam_scan.txt max_e_balue min_non_N_CDS_len";

$gff3_file = $ARGV[0];
$fasta_file = $ARGV[1];
$pfam_file = $ARGV[2];
$max_e_value = $ARGV[3];
$min_CDS_len = $ARGV[4];


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

open($in, $fasta_file);
while (($name, $s) = fasta_nonwhite_get($in)) {
	$seq{$name} = $s;
}
close $in;

while (($mRNA_id, $info) = each %gene) {
	@pos = sort {$a <=> $b} @{$info->{pos}};
	$CDS = '';
	for ($i = 0; $i < $#pos; $i += 2) {
		$CDS .= substr($seq{$info->{name}}, $pos[$i] - 1, $pos[$i + 1] - $pos[$i] + 1);
	}

	if ($info->{strand} eq '-') {
		$CDS = rev_comp($CDS);
	}

	$non_N_len = ($CDS =~ tr/ATGCatgc/ATGCatgc/);
	if ($non_N_len >= $min_CDS_len) {
		$gene_id = gff3_tree_get_root_id(\%gff3_tree, $mRNA_id);
		$len_pass_flag{$gene_id} = 1;
	}
}

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

	if ($len_pass_flag{$gene_id} or $pfam_pass_flag{$gene_id}) {
		print "$_\n";
	}
}
close $in;
