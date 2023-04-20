#!/usr/bin/perl

use my_bio;
use my_gff;

(@ARGV != 2) and die "usage: $0 gff3 fasta\n";

gff3_tree_set_relations(\%gff3_tree, $ARGV[0]);

open IN, $ARGV[0];
while (<IN>) {
	chomp;
	($seq_name, $feature, $st, $ed, $strand, $attr) = (split(/\t/, $_))[0, 2, 3, 4, 6, 8];
	$id = gff3_attr_get_id($attr);

	if ($feature eq 'exon') {
		@gene_ids = gff3_tree_get_parent_id(\%gff3_tree, $id);
		for $gene_id (@gene_ids) {
			if (not $gene_order_flag{$gene_id}) {
				push(@gene_order, $gene_id);
				$gene_order_flag{$gene_id} = 1;
			}
			$gene{$gene_id}{name} = $seq_name;
			$gene{$gene_id}{strand} = $strand;
			push(@{$gene{$gene_id}{pos}}, [$st, $ed]);
		}
	}
	elsif ($feature eq 'mRNA') {
		if ($attr =~ /Note=([^;]+)/) {
			$gene_note{$id} = $1;
			$gene_note{$id} =~ s/,/ /g;
			$gene_note{$id} =~ s/%(..)(?{ chr(hex($1)) })/$^R/g;
		}
	}
}
close IN;

open($in, $ARGV[1]);
while (($name, $s) = fasta_nonwhite_get($in)) {
	$seq{$name} = $s;
}
close $in;


for $gene_id (@gene_order) {
	$info = $gene{$gene_id};
	@pos = sort {$a->[0] <=> $b->[1]} @{$info->{pos}};
	$mRNA = '';
	for $st_ed (@pos) {
		$mRNA .= substr($seq{$info->{name}}, $st_ed->[0] - 1, $st_ed->[1] - $st_ed->[0] + 1);
	}

	if ($info->{strand} eq '-') {
		$mRNA = rev_comp($mRNA);
	}

	if (defined $gene_note{$gene_id}) {
		printf(">%s %s\n", $gene_id, $gene_note{$gene_id});
	}
	else {
		printf(">%s\n", $gene_id);
	}
	print_seq($mRNA);
}
