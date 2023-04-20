#!/usr/bin/perl

use strict;
use warnings;
use constant WIDTH => 80;

&main;


sub main
{
    my ( %gene, %cds, %exon, %fasta, %tree );
    %fasta = &input_fasta;
    &set_tree( \%tree );
    &input_gff( \%gene, \%cds, \%exon, \%tree );
    &cal_stats( \%gene, \%cds, \%exon, \%fasta )
}


sub set_tree
{
    my ( $tree ) = @_;
    my @buf;
    my ( $ID, $parent );
    open IN, $ARGV[0];
    while( <IN> )
    {
	chomp;
	next if( $_ =~ /^#/ );
	next if( $_ eq "" );
	@buf = split/\t/;
	if( $buf[2] eq "gene" )
	{
	    next if( $buf[8] !~ /ID=(\S+)/ );
	    $tree->{( split( /;/, $1 ) )[0]} = "";
	}else{
	    next if( $buf[8] !~ /ID=(\S+);Parent=(\S+)/ );
	    $parent = (split( /;/, $2 ))[0];
	    $tree->{$1} = $parent;
	}
    }
    close IN;
}

sub input_fasta
{
    my %fasta;
    my $ID;
    open IN, $ARGV[1];
    while( <IN> )
    {
	chomp;
	if( $_ =~ />(\S+)/ )
	{
	    $ID = $1;
	}else{
	    $fasta{$ID} .= $_;
	}
    }
    close IN;
    return %fasta;
}


sub input_gff
{
    my ( $gene, $cds, $exon, $tree ) = @_;
    my @buf;
    my ( $ID, $parent, $gene_ID );
    open IN, $ARGV[0];
    while( <IN> )
    {
	chomp;
	next if( $_ =~ /^#/ );
	next if( $_ eq "" );
	@buf = split/\t/;
	if( $buf[2] eq "gene" ){
	    next if( $buf[8] !~ /ID=(\S+)/ );
	    $ID = ( split( /;/, $1 ) )[0];
	    $gene->{ $ID } = [ @buf ];
	}elsif( $buf[2] eq "exon" ){	   
	    $buf[8] =~ /ID=(\S+);Parent/;
	    ( $gene_ID, $parent ) = &set_parent( \$1, \%{$tree} );
	    #print $buf[2],"\t".$1,"\t".$parent,"\n";
	    push( @{$exon->{ $gene_ID }{ $parent }}, [ @buf ] );
	}elsif( $buf[2] eq "CDS" ){
	    $buf[8] =~ /ID=(\S+);Parent/;
	    ( $gene_ID, $parent ) = &set_parent( \$1, \%{$tree} );
	    #print $buf[2],"\t".$1,"\t".$parent,"\n";
	    push( @{$cds->{ $gene_ID }{ $parent }}, [ @buf ] );
	}
    }
    close IN;
}

sub set_parent
{
    my ( $ID, $tree ) = @_;
    my $parent = $tree->{$$ID};
    my $gene_ID = $parent;
    while( 1 )
    {
	last if( ! exists $tree->{$gene_ID} );
	last if( $tree->{$gene_ID} eq "" );
	$gene_ID = $tree->{$gene_ID};
    }
    return ( $gene_ID, $parent );
}


sub cal_stats
{
    my ( $gene, $cds, $exon, $seq ) = @_;
    
    my $gene_num = 0;
    my $total_ex_int_len = 0;    

    my $single_exon = 0;
    my $total_exon_len = 0;
    my $total_exon_num = 0;

    my $total_CDS_len = 0;
    my $total_CDS_num = 0;

    my $total_intron_num = 0;
    my $total_intron_len = 0;
    my $total_GTAG_site = 0;
    my $splice_site;
    my %splice = ( "GTAG" => 1, "CTAC" => 1 );
    foreach my $gene ( sort keys %{$cds} )
    {
	my %len = ();
	next if( ! exists $exon->{$gene} );
	foreach my $prod ( sort keys %{$cds->{$gene}} )
	{
	    next if( ! exists $exon->{$gene}{$prod} );
	    next if( ! exists $exon->{$gene}{$prod} );
	    for( my $i = 0; $i <= $#{$cds->{$gene}{$prod}}; $i++ )
	    {
		my @buf = @{$cds->{$gene}{$prod}[$i]};
		$len{$prod}[0] += $buf[4] - $buf[3] + 1;
	    }
	    for( my $i = 0; $i <= $#{$exon->{$gene}{$prod}}; $i++ )
	    {
		my @buf = @{$exon->{$gene}{$prod}[$i]};
		$len{$prod}[2] += $buf[4] - $buf[3] + 1;
		$len{$prod}[1]++;
	    }
	}
	foreach my $longest ( sort { $len{$b}[0] <=> $len{$a}[0] or $len{$a}[1] <=> $len{$b}[1] or $a cmp $b } keys %len )
	{
	    $gene_num++;
	    $total_CDS_len += $len{$longest}[0];
	    $total_exon_num += $len{$longest}[1];
	    $total_exon_len += $len{$longest}[2];
	    $single_exon++ if( $len{$longest}[1] == 1 );
	    my @buf = sort { $a->[3] <=> $b->[3] } @{$exon->{$gene}{$longest}};
	    $total_ex_int_len += $buf[$#buf][4] - $buf[0][3] + 1;
	    $total_intron_num += scalar @buf - 1;
	    for( my $i = 0; $i <= $#buf-1; $i++ )
	    {
		next if( $buf[$i][0] ne $buf[$i+1][0] );
		$total_intron_len += $buf[$i+1][3] - 1 - ( $buf[$i][4] + 1 ) + 1;
		$splice_site = substr( $seq->{$buf[$i][0]}, $buf[$i][4], 2 ).substr( $seq->{$buf[$i+1][0]}, $buf[$i+1][3] - 3, 2 );
		$splice_site =~ tr/a-z/A-Z/;
		$total_GTAG_site++ if( exists $splice{ $splice_site } );
	    }
	    #print $gene,"\t".$longest,"\t".$len{$longest}[1],"\t".$len{$longest}[0],"\t".$len{$longest}[2],"\n";
	    last;
	}
    }
    print "#gene:\t".$gene_num,"\n";
    print "total exon+intron length:\t".$total_ex_int_len,"\n";
    print "mean exon+intron length:\t".sprintf( "%.2f", $total_ex_int_len/$gene_num ),"\n\n"; 
    
    print "#exons:\t".$total_exon_num,"\n";
    print "#exons/#genes:\t".sprintf( "%.2f", $total_exon_num/$gene_num ),"\n";
    print "#single exon genes:\t".$single_exon,"\n";
    print "total exon length:\t".$total_exon_len,"\n";
    print "mean exon length:\t".sprintf( "%.2f", $total_exon_len /$total_exon_num ),"\n";
    print "mean CDS length:\t".sprintf( "%.2f", $total_CDS_len / $gene_num ),"\n\n";

    print "#intron:\t".$total_intron_num,"\n";
    print "total intron length:\t".$total_intron_len,"\n";
    if( $total_intron_num != 0 ){
	print "mean intron length:\t".sprintf( "%.2f", $total_intron_len / $total_intron_num ),"\n";
	print "GT-AG splicing site:\t".sprintf( "%.2f", $total_GTAG_site * 100/ $total_intron_num ),"\n";
    }else{
	print "mean intron length:\t---\n";
	print "GT-AG splicing site:\t---\n";
    }

}
