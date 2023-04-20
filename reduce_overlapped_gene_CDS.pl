#!/usr/bin/perl

use my_bio;
use List::Util qw/min/;

(@ARGV != 3) and die "usage: $0 genes.gff3 ref.fasta min_overlap_rate\n";

$min_ol_rate = $ARGV[2];

for $i (0..1) {
	gff3_tree_set_relations(\%{$gff3_tree[$i]}, $ARGV[0]);
}

set_fasta_seq(\%ref_seq, $ARGV[1]);
for $i (0..1) {
	set_gene_info(\%{$gene_info[$i]}, $ARGV[0], \%{$gff3_tree[$i]}, \%ref_seq);
}

set_replace_list($ARGV[0], $ARGV[0], \%{$gene_info[0]}, \%{$gene_info[1]}, \%replace_list, \%{$gff3_tree[0]}, \%{$gff3_tree[1]}, $min_ol_rate);

print_unreplaced_gene($ARGV[0], \%replace_list, \%{$gff3_tree[0]});


sub set_fasta_seq
{
	my $ref_seq = shift;
	my $file_name = shift;
	my ($in, $name, $seq);

	open($in, $file_name);
	while (($name, $seq) = fasta_nonwhite_get($in)) {
		$ref_seq->{$name} = $seq;
	}
	close $in;
}


sub set_gene_info 
{
	my $gene_info = shift;
	my $file = shift;
	my $gff3_tree = shift;
	my $ref_seq = shift;

	my($in, $l, $feat, $st, $ed, $attr, $mRNA_id, $gene_id, $seq_name);

	open($in, $file);
	while (chomp($l = <$in>)) {
		($seq_name, $feat, $st, $ed, $strand, $attr) = (split(/\t/, $l))[0, 2, 3, 4, 6, 8];

		if ($feat eq 'CDS') {
			$gene_id  = gff3_tree_get_root_id($gff3_tree, gff3_attr_get_id($attr));
			$gene_info->{$gene_id}[0] += (substr($ref_seq->{$seq_name}, $st - 1, $ed - $st + 1) =~ tr/ATGCatgc/ATGCatgc/);
			$gene_info->{$gene_id}[1] = $strand;
		}
		elsif ($feat eq 'exon') {
			$gene_id = gff3_tree_get_root_id($gff3_tree, gff3_attr_get_id($attr));
			$gene_info->{$gene_id}[2] += ($ed - $st + 1);
		}

	}
	close $in;
}

sub set_replace_list
{
	my $file1 = shift;
	my $file2 = shift;
	my $gene_info1 = shift;
	my $gene_info2 = shift;
	my $replace_list = shift;
	my $gff3_tree1 = shift;
	my $gff3_tree2 = shift;
	my $min_ol_rate = shift;

	my($in, $l, $feat1, $attr1, $feat2, $attr2, $gene_id1, $gene_id2, %exon_ol_len, $ol_len, $key);

	open($in, "bedtools intersect -a $file1 -b $file2 -wao |");
	while (chomp($l = <$in>)) {
		($feat1, $st1, $ed1, $attr1, $feat2, $st2, $ed2,  $attr2, $ol_len) = (split(/\t/, $l))[2, 3, 4, 8, 11, 12, 13, 17, 18];
#		next unless ($feat1 eq 'exon' and $feat2 eq 'exon');
		next unless ($feat1 eq 'CDS' and $feat2 eq 'CDS');

		$gene_id1  = gff3_tree_get_root_id($gff3_tree1, gff3_attr_get_id($attr1));
		$gene_id2  = gff3_tree_get_root_id($gff3_tree2, gff3_attr_get_id($attr2));

		if ($gene_id1 ne $gene_id2) {
			$exon_ol_len{$gene_id1 . ';' . $gene_id2} += $ol_len;
		}
	}
	close $in;

	while (($key, $ol_len) = each %exon_ol_len) {
		($gene_id1, $gene_id2) = split(/;/, $key);

#		if ($gene_info1->{$gene_id1}[1] eq $gene_info2->{$gene_id2}[1] and $ol_len / min(($gene_info1->{$gene_id1}[2], $gene_info2->{$gene_id2}[2])) >= $min_ol_rate) {
		if ($gene_info1->{$gene_id1}[1] eq $gene_info2->{$gene_id2}[1] and $ol_len / min(($gene_info1->{$gene_id1}[0], $gene_info2->{$gene_id2}[0])) >= $min_ol_rate) {
			if ($gene_info1->{$gene_id1}[0] < $gene_info2->{$gene_id2}[0]) {
				$replace_list->{$gene_id1} = $gene_id2;
			}
			elsif ($gene_info1->{$gene_id1}[0] > $gene_info2->{$gene_id2}[0]) {
				$replace_list->{$gene_id2} = $gene_id1;
			}
			else {
				if ($gene_info1->{$gene_id1}[2] < $gene_info2->{$gene_id2}[2]) {
					$replace_list->{$gene_id1} = $gene_id2;
				}
				elsif ($gene_info1->{$gene_id1}[2] > $gene_info2->{$gene_id2}[2]) {
					$replace_list->{$gene_id2} = $gene_id1;
				}
				else {
					if (not($replace_list->{$gene_id1}) and not($replace_list->{$gene_id2})) {
						$replace_list->{$gene_id2} = $gene_id1;
					}
				}
			}
		}
	}
}

sub print_unreplaced_gene
{
	my $file = shift;
	my $replace_list = shift;
	my $gff3_tree = shift;

	my($in, $l, $feat, $st, $ed, $attr, $gene_id);

	open($in, $file);
	while (chomp($l = <$in>)) {
		if ($l eq '' or $l =~ /^#/) {
			next;
		}
		($feat, $st, $ed, $attr) = (split(/\t/, $l))[2, 3, 4, 8];

		$gene_id  = gff3_tree_get_root_id($gff3_tree, gff3_attr_get_id($attr));
		
		unless ($replace_list->{$gene_id}) {
			print "$l\n";
		}
	}
	close $in;
}

sub gff3_tree_set_relations
{
	my ($in, $source, $feat, $st, $ed, $attr, $id, $parent, %tree);

	my $tree_ref = shift;
	my $file_name = shift;

	open($in, $file_name);
	while (chomp($l = <$in>)) {
		($source, $feat, $st, $ed, $attr) = (split(/\t/, $l))[1, 2, 3, 4, 8];
		if ($attr =~ /ID=([^;]+)/) {
			$id = $1;
		}
		else {
			next;
		}

		if ($attr =~ /Parent=([^;]+)/) {
			$parent = $1;
		}
		else {
			$parent = '';
		}

		$tree_ref->{$id} = [split(/,/, $parent)];
	}
	close $in;
}

sub gff3_tree_get_parent_id
{
	my $tree_ref = shift;
	my $id = shift;

	return @{$tree_ref->{$id}};
}

sub gff3_tree_get_root_id
{
	my $tree_ref = shift;
	my $id = shift;

	while ($tree_ref->{$id}[0]) {
		$id = $tree_ref->{$id}[0];
	}

	return $id;
}

sub gff3_attr_get_id
{
	my $attr = shift;
	my $id;
	if ($attr =~ /ID=([^;]+)/) {
		$id = $1;
	}
	return $id;
}
