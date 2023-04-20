#!/usr/bin/perl

(@ARGV != 2) and die "usage: $0 gff3 mRNA_id_desc.tsv\n";

gff3_tree_set_relations(\%gff3_tree, $ARGV[0]);

open($in, $ARGV[1]);
while (chomp($l = <$in>)) {
	@f = split(/\t/, $l);
	$mRNA_id2desc{$f[0]} = [@f[1..$#f]]
}
close $in;

open($in, $ARGV[0]);
while (chomp($l = <$in>)) {
	($feat, $attr) = (split(/\t/, $l))[2, 8];
	$id = gff3_attr_get_id($attr);

	if ($feat eq 'mRNA' and defined $mRNA_id2desc{$id}) {
		$gene_id = gff3_tree_get_root_id(\%gff3_tree, $id);
		push(@{$gene_id2desc{$gene_id}}, @{$mRNA_id2desc{$id}});
	}
}
close $in;

open($in, $ARGV[0]);
while (chomp($l = <$in>)) {
	($feat, $attr) = (split(/\t/, $l))[2, 8];
	$id = gff3_attr_get_id($attr);

	if ($feat eq 'mRNA') {
		if (defined $mRNA_id2desc{$id}) {
			if ($l !~ /;$/) {
				$l .= ';';
			}
			$l .= sprintf("Note=%s", join(',', map({gff3_excape_attribute_character($_)} @{$mRNA_id2desc{$id}})));
		}
	}
	elsif ($feat eq 'gene') {
		if (defined $gene_id2desc{$id}) {
			if ($l !~ /;$/) {
				$l .= ';';
			}
			$l .= sprintf("Note=%s", join(',', map({gff3_excape_attribute_character($_)} @{$gene_id2desc{$id}})));
	
		}
	}

	print "$l\n";
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


BEGIN {
	my %escape = (
		',' => '%2C',
		';' => '%3B',
		'=' => '%3D',
	);

	sub gff3_excape_attribute_character
	{
		my $attr = shift;	
		
		for ($i = length($attr) - 1; $i >= 0; --$i) {
#print(join("\t", (substr($attr, $i, 1), $escape{substr($attr, $i, 1)})), "\n");
			if (defined $escape{substr($attr, $i, 1)}) { 
				substr($attr, $i, 1, $escape{substr($attr, $i, 1)});
			}
		}

		return $attr;
	}
}
