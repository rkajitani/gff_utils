#!/usr/bin/perl

use my_bio;
use my_gff;

(@ARGV != 3) and die "usage: $0 genes.gff3 genome.fa min_non_N_CDS_len";

$min_CDS_len = $ARGV[2];

gff3_tree_set_relations(\%gff3_tree, $ARGV[0]);

open($in, $ARGV[0]);
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

open($in, $ARGV[1]);
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
		$pass_flag{$gene_id} = 1;
	}
}

open($in, $ARGV[0]);
while (<$in>) {
	chomp;
	($attr) = (split(/\t/, $_))[8];

	$id = gff3_attr_get_id($attr);
	$gene_id = gff3_tree_get_root_id(\%gff3_tree, $id);

	if ($pass_flag{$gene_id}) {
		print "$_\n";
	}
}
close $in;
