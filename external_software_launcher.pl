#!/usr/bin/perl -w

use strict;
use warnings;

my $type       = shift;
my $low_bound  = shift;
my $high_bound = shift;
my $step       = shift;
my $directory  = shift; #where to find the new family sequences (.fasta files)
my $suffix     = shift;
my $output_dir = shift;
my $ref_hmm_dir = shift;
print "Launching $type $low_bound $high_bound $step\n";

my $tree_dir = $output_dir."/trees";
my $hmm_dir  = $output_dir."/HMMs";
my $aln_dir  = $output_dir."/alignments";

##if the $output_dir does not exists create it.
print "Creating $output_dir\n" unless -e $output_dir;
`mkdir -p $output_dir`         unless -e $output_dir;
##if the $tree_dir does not exists create it.
print "Creating $tree_dir\n" unless -e $tree_dir;
`mkdir -p $tree_dir`         unless -e $tree_dir;
##if the $hmm_dir does not exists create it.
print "Creating $hmm_dir\n" unless -e $hmm_dir;
`mkdir -p $hmm_dir`         unless -e $hmm_dir;
##if the $aln_dir does not exists create it.
print "Creating $aln_dir\n" unless -e $aln_dir;
`mkdir -p $aln_dir`         unless -e $aln_dir;

if ( $type eq 'new' ) {
	for ( my $i = $low_bound; $i <= $high_bound; $i = $i + $step ) {
		print "processing family : $i\t$directory/$i$suffix\n";
		
		if ( -e "$directory/$i$suffix" ) {
			##align the fasta file
			`muscle -in $directory/$i$suffix -out $aln_dir/$i.aln` unless -e "$directory/$i.aln";
			##HMM of the alignment file
			if (-e "$hmm_dir/$i.hmm.gz" && !-e "$hmm_dir/$i.hmm"){
				##gunzip the file
				`gunzip $hmm_dir/$i.hmm.gz`;
			}
			`hmmbuild --informat afa $hmm_dir/$i.hmm $aln_dir/$i.aln` unless -e "$hmm_dir/$i.hmm";
			if (-e "$hmm_dir/$i.hmm" && !-e "$hmm_dir/$i.hmm.gz"){
				##zip the file
				`gzip $hmm_dir/$i.hmm`;
			}
			##Tree of the alignment
			`fasttree $aln_dir/$i.aln > $tree_dir/$i.tree` unless -e "$tree_dir/$i.tree";
		}
	}
} elsif ( $type eq 'old' ) {
	for ( my $i = $low_bound; $i <= $high_bound; $i = $i + $step ) {
		##align sequences
		`hmmalign --outformat afa --trim --amino -o $aln_dir/$i.aln $ref_hmm_dir/$i.aln $directory/$i.aln`;
		##HMM of the alignment file
		`hmmbuild --informat afa $hmm_dir/$i.hmm $aln_dir/$i.aln` unless -e "$hmm_dir/$i.hmm";
		##Tree of the alignment
		`fasttree $aln_dir/$i.aln > $tree_dir/$i.tree` unless -e "$tree_dir/$i.tree";
	}
} else {
	die "Command $type not recognized\n";
}
