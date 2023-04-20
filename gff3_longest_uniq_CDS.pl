#!/usr/bin/perl

(@ARGV != 1) and die "usage: $0 gff3\n";


gff3_tree_set_relations(\%gff3_tree, $ARGV[0]);

open($in, $ARGV[0]);
while (chomp($l = <$in>)) {
	($source, $feat, $st, $ed, $attr) = (split(/\t/, $l))[1, 2, 3, 4, 8];
	$id = gff3_attr_get_id($attr);

	if ($feat eq 'CDS') {
		for $mRNA_id (gff3_tree_get_parent_id(\%gff3_tree, $id)) {
			$CDS_len{$mRNA_id} += ($ed - $st);
		}
	}
	elsif ($feat eq 'gene') {
		$gene_pos{$id} = $st . '_' . $ed;
		$gene_desc_line{$id} = "$l\n";
	}
}
close $in;

for $mRNA_id (keys %CDS_len) {
	$pos = $gene_pos{gff3_tree_get_root_id(\%gff3_tree, $mRNA_id)};

	if ($CDS_len{$mRNA_id} > $max_CDS_len{$pos}) {
		$max_CDS_len{$pos} = $CDS_len{$mRNA_id};
		$max_CDS_mRNA_id{$pos} = $mRNA_id;
	}
}

open($in, $ARGV[0]);
while (chomp($l = <$in>)) {
	($source, $feat, $attr) = (split(/\t/, $l))[1, 2, 8];

	if ($feat eq 'gene') {
		next;
	}
	elsif ($feat eq 'mRNA') {
		@mRNA_ids = (gff3_attr_get_id($attr));
	}
	else {
		$id = gff3_attr_get_id($attr);
		@mRNA_ids = gff3_tree_get_parent_id(\%gff3_tree, $id);
	}
	$gene_id = gff3_tree_get_root_id(\%gff3_tree, $mRNA_ids[0]);

	for $mRNA_id (@mRNA_ids) {
		if ($mRNA_id eq $max_CDS_mRNA_id{$gene_pos{$gene_id}}) {
			print $gene_desc_line{$gene_id};
			$gene_desc_line{$gene_id} = '';

			if (@mRNA_ids > 1) {
				@f = split(/\t/, $l);
				$f[8] =~ s/Parent=([^;]+)/Parent=$mRNA_id/;
				$l = join("\t", @f);
			}
			print "$l\n";
		}
	}
}
close $in;


sub gff3_tree_set_relations
{
	my ($in, $source, $feat, $st, $ed, $attr, $id, $parent, %tree);

	my $tree_ref = shift;
	my $file_name = shift;

	open($in, $ARGV[0]);
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
