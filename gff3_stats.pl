#!/usr/bin/perl

(@ARGV != 1) and die "usage: $0 gene.gff3\n";

#$out_prefix = $ARGV[1];

$bin = 100;

open IN, $ARGV[0];
while (<IN>) {
	($seq_name, $feature, $st, $ed, $attr) = (split(/\t/, $_))[0, 2, 3, 4, 8];
	next if ($feature ne 'exon' or $attr =~ /tRNA/);

	$attr =~ /Parent=([^;]+)/;
	$gene_id = $1;
	$gene_id =~ s/\|m.*$//;
	$gene_id =~ s/-mRNA-\d+$//;

	++$stats{$gene_id}{num_exon};
	$stats{$gene_id}{len_exon} += $ed - $st + 1;
	$stats{$gene_id}{st} = $st if (not defined($stats{$gene_id}{st}) or $st <= $stats{$gene_id}{st});
	$stats{$gene_id}{ed} = $ed if ($ed > $stats{$gene_id}{ed});
	++$exon_len_hist[int(($ed - $st + 1) / $bin)];
}
close IN;


for (values %stats) {
	++$num_gene;
	$num_exon += $_->{num_exon};
	$total_len_exon += $_->{len_exon};
	$total_len_gene += $_->{ed} - $_->{st} + 1;
	++$num_single_exon if ($_->{num_exon} == 1);
	++$gene_len_hist[int(($_->{ed} - $_->{st} + 1) / $bin)];
}

#print "num_gene = $num_gene\n";
#print "num_exon = $num_exon\n";
#print("num_exon/gene = ", $num_exon / $num_gene, "\n");
#print "num_single_exon = $num_single_exon\n";
#print "total_len_gene = $total_len_gene\n";
#print "total_len_exon = $total_len_exon\n";
#print("mean_len_gene = ", $total_len_gene / $num_gene, "\n");
#print("mean_len_exon = ", $total_len_exon / $num_exon, "\n");

print "$num_gene\n";
print "$total_len_gene\n";
print($total_len_gene / $num_gene, "\n");
print "$num_exon\n";
print($num_exon / $num_gene, "\n");
print "$num_single_exon\n";
print "$total_len_exon\n";
print($total_len_exon / $num_exon, "\n");


#print "\n#gene length histogram\n";
#for (0..$#gene_len_hist) {
#	print(($_ * $bin + $bin/2), "\t", $gene_len_hist[$_] + 0, "\n");
#}
#
#print "\n#exon length histogram\n";
#for (0..$#exon_len_hist) {
#	print(($_ * $bin + $bin/2), "\t", $exon_len_hist[$_] + 0, "\n");
#}
