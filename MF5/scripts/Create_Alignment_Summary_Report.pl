#!/usr/bin/perl
# Priscila Darakjian
# 
# This script extracts specific information from 
# each result summary file from each sample for which 
# alignment was run through STAR
#
# Usage: perl Create_Report.pl pathtoresultslocdir
#

my $results_dir = @ARGV[0]; 
opendir(DIR,"$results_dir") or die "$!";
open(OUT,">STAR_Results_Summary.csv") or die "$!";

my @files = grep {/Log\.final\.out$/} readdir DIR;
close DIR;
my $first_file=1;

foreach my $file(@files){
	open(IN,"$results_dir"."/$file") or die "$!";
	print "$file\n";
	if ($first_file==1){
		print OUT "Sample,Started,Started mapping,Finished,Mapping speed ( Million of reads per hour),,Total Reads,Avg Length,UNIQUE READS:, Uniquely mapped,Uniquely mapped %,Average mapped length,Number of splices: Total,Number of splices: Annotated (sjdb),Number of splices: GT/AG,Number of splices: GC/AG,Number of splices: AT/AC,Number of splices: Non-canonical,Mismatch rate per base (%),Deletion rate per base,Deletion average length,Insertion rate per base,Insertion average length,MULTI-MAPPING READS:,Number of reads mapped to multiple loci,% of reads mapped to multiple loci,Number of reads mapped to too many loci,% of reads mapped to too many loci,UNMAPPED READS:,% of reads unmapped: too many mismatches,% of reads unmapped: too short,% of reads unmapped: other";
		$first_file=0;
	}
	my $sample = substr($file,0,-13);
	print OUT "\n" . $sample;
	while(<IN>){
		my $row = $_;
		my @line=split('\t',$row);
		chomp(@line);
		#my $item = $line[0];
		my $value = $line[1];
		print OUT "," . $line[1];
	}
}
