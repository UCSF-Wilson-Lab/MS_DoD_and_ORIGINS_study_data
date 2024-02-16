#!/usr/bin/perl
use strict;
use warnings 'FATAL' => 'all';

# PATH to our local project directory is not provided
my $project_dir = "[PATH to project folder]"
my $results_dir = $project_dir . "/resources/fig2_CDHIT_results/";
my $cdhit_results_fh = $results_dir . "cluster_output_motif_and_alt_species.clstr";
# OUTPUT
my $results_fh       = $results_dir . "cluster_output_motif_and_alt_species.FORMATTED.clstr";

my %clusters;  # has to store all sequence cluster info
my %sizes;     # number of peptides per cluster
my $current_cluster;


### 1. Store Cluster info
open(my $IN, "<", $cdhit_results_fh) or die "ERROR: Cannot read file\n";
while (my $line = <$IN>){
	chomp($line);
	if($line =~ /^>/){
		my $clusid = (split /\s+/,$line)[-1];
		$current_cluster = $clusid;
	} else{
		my $iter    = (split /\s+/,$line)[0];
		my $peptide = (split />/,$line)[1];
		$peptide = (split /\.\.\./,$peptide)[0];
		my $seqident = (split /\s+/,$line)[-1];
		$seqident =~ s/%//;
		my $entry = "$current_cluster\t$peptide\t$seqident";
		$clusters{$current_cluster}{$iter} = $entry;
		$sizes{$current_cluster}++;
	}
	
}
close $IN;


### 2. OUTPUT formatted table
open(my $OUT, ">",$results_fh) or die "ERROR: Cannot write file\n";
print $OUT "ClusterID\tPeptide\t%Ident\tClusterSize\n";

foreach my $clus (keys %clusters){
	foreach my $iter (keys %{$clusters{$clus}}){
		my $entry     = $clusters{$clus}{$iter};
		my $clus_size = $sizes{$clus};
		
		$entry .= "\t$clus_size";
		print $OUT "$entry\n";
	}
}
close $OUT;

print "\n\n>>> DONE <<<\n\n";