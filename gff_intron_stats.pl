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
	if ($feature eq 'exon' and $attr =~ /Parent=([^;]*)/) {
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
	if (@{$_->{pos}} > 2) {
		@pos = sort{$a <=> $b} @{$_->{pos}};
		for ($i = 1; $i < @pos - 2; $i += 2) {
			++$num_intron;
#			$intron_len_hist[int(($pos[$i + 1] - $pos[$i] + 1) / $bin)];
#			$total_len_intron += ($pos[$i + 1] - $pos[$i] + 1);
			$intron_len_hist[int(($pos[$i + 1] - $pos[$i] - 1) / $bin)];
			$total_len_intron += ($pos[$i + 1] - $pos[$i] - 1);

			$s = substr($$seq, $pos[$i], 2) . substr($$seq, $pos[$i + 1] - 3, 2); 
			if ($_->{strand} eq '-') {
				$s = rev_comp($s);
			}
			++$splice_site{$s};
		}
	}
}

print "num_intron = $num_intron\n";
print "total_len_intron = $total_len_intron\n";
print("mean_len_intron = ", $total_len_intron / $num_intron, "\n");
print("GTAG_splice_site = ", $splice_site{GTAG} / $num_intron * 100, "\n");

print "$num_intron\n";
print "$total_len_intron\n";
print($total_len_intron / $num_intron, "\n");
print($splice_site{GTAG} / $num_intron * 100, "\n");


#print "\n#gene length histogram\n";
#for (0..$#gene_len_hist) {
#	print(($_ * $bin + $bin/2), "\t", $gene_len_hist[$_] + 0, "\n");
#}
#
#print "\n#intron length histogram\n";
#for (0..$#intron_len_hist) {
#	print(($_ * $bin + $bin/2), "\t", $intron_len_hist[$_] + 0, "\n");
#}
