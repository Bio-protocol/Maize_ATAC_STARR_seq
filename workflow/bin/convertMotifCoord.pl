#!/usr/bin/perl
use strict;
use warnings;

# iterate over GFF file
open F, $ARGV[0] or die;
while(<F>){
	chomp;

	# skip header
	if($_ =~ /^##/){
		next;
	}

	# split row by tab
	my @col = split("\t",$_);

	# correct position
	my @coord = split(":",$col[0]);
	my $chr = $coord[0];
	my @pos = split("-",$coord[1]);
	my $str = $pos[0];
	my $motif_str = $col[3] - 1;
	my $motif_end = $col[4];
	my $adjstr = $str + $motif_str;
	my $adjend = $str + $motif_end;

	# get p-value
	my @fields = split(";",$col[8]);
	my $pval = 1;
	my $motifID;
	foreach(@fields){
		if($_ =~ /^Name/){
			my @id = split("=",$_);
			my @idd = split("_",$id[1]);
			$motifID = join("_",$idd[0],$idd[1]);
		}elsif($_ =~ /^pvalue/){
			my @vals = split("=",$_);
			$pval = $vals[1];
		}
	}

	# print output
	if($pval < 0.00001){
		print "$chr\t$adjstr\t$adjend\t$motifID\n";
	}
}
close F;
