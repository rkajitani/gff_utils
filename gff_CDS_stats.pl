#!/usr/bin/perl

use my_bio;

(@ARGV != 2) and die "usage: $0 gene.gff3 genome.fa\n";

open($in, $ARGV[1]);
while (($name, $seq) = fasta_nonwhite_get($in)) {
	$ref{$name} = $seq;
}
close $in;

$bin = 100;

open IN, $ARGV[0];
while (<IN>) {
	chomp;
	($seq_name, $feature, $st, $ed, $strand, $attr) = (split(/\t/, $_))[0, 2, 3, 4, 6, 8];
	if ($feature eq 'CDS' and $attr =~ /ID=([^;]*)/) {
		$id = $1;
		$mRNA{$id}{seq_name} = $seq_name;
		$mRNA{$id}{strand} = $strand;
		push(@{$mRNA{$id}{pos}}, $st);
		push(@{$mRNA{$id}{pos}}, $ed);
	}
}
close IN;


for (values %mRNA) {
	$seq = \$ref{$_->{seq_name}};
	if (@{$_->{pos}} >= 2) {
		@pos = sort{$a <=> $b} @{$_->{pos}};
		$s = '';
		for ($i = 0; $i < @pos - 1; $i += 2) {

			$s .= substr($$seq, $pos[$i] - 1, ($pos[$i + 1] - $pos[$i] + 1)); 
		}
		if ($_->{strand} eq '-') {
			$s = rev_comp($s);
		}
		++$num_CDS;
		$CDS_len_hist[int(length($s) / $bin)];
		$total_len_CDS += length($s);
	}
}

print "num_CDS = $num_CDS\n";
#print "total_len_CDS = $total_len_CDS\n";
print("mean_len_CDS = ", $total_len_CDS / $num_CDS, "\n");

print "$num_CDS\n";
#print "$total_len_CDS\n";
print($total_len_CDS / $num_CDS, "\n");


#print "\n#gene length histogram\n";
#for (0..$#gene_len_hist) {
#	print(($_ * $bin + $bin/2), "\t", $gene_len_hist[$_] + 0, "\n");
#}
#
#print "\n#CDS length histogram\n";
#for (0..$#CDS_len_hist) {
#	print(($_ * $bin + $bin/2), "\t", $CDS_len_hist[$_] + 0, "\n");
#}
