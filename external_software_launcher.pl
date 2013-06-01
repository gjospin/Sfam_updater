#!/usr/bin/perl -w
use strict;
use warnings;
use Bio::AlignIO;

my $type        = shift;
my $low_bound   = shift;
my $high_bound  = shift;
my $step        = shift;
my $directory   = shift;    #where to find the new family sequences (.fasta files)
my $suffix      = shift;
my $output_dir  = shift;
my $ref_hmm_dir = shift;

my %fc_range=();
open(FCIN,"$ref_hmm_dir/fc_range.txt");
while(<FCIN>){
	chomp($_);
	next if $_ =~ m/^#/;
	my @line = split(/\t/,$_);
	# line should be formatted as FC_ID<tab>low_family_id<tab>high_family_id
	push(@{$fc_range{$line[0]}} ,($line[1],$line[2]));
}
close(FCIN);

print "Launching $type $low_bound $high_bound $step\n";
my $tree_dir      = $output_dir."/trees";
my $hmm_dir       = $output_dir."/HMMs";
my $aln_dir       = $output_dir."/alignments";
my $rep_dir       = $directory."/representatives";
my $seed_hmm_dir  = $output_dir."/seedHMMs";
my $aln_stock_dir = $output_dir."/alignments_stock";
##if the $output_dir does not exists create it.
#print "Creating $output_dir\n" unless -e $output_dir;
#`mkdir -p $output_dir`         unless -e $output_dir;
##if the $tree_dir does not exists create it.
#print "Creating $tree_dir\n" unless -e $tree_dir;
#`mkdir -p $tree_dir`         unless -e $tree_dir;
##if the $hmm_dir does not exists create it.
#print "Creating $hmm_dir\n" unless -e $hmm_dir;
#`mkdir -p $hmm_dir`         unless -e $hmm_dir;
##if the $aln_dir does not exists create it.
#print "Creating $aln_dir\n" unless -e $aln_dir;
#`mkdir -p $aln_dir`         unless -e $aln_dir;
##if the $rep_dir does not exists create it.
#print "Creating $rep_dir\n" unless -e $rep_dir;
#`mkdir -p $rep_dir`         unless -e $rep_dir;
##if the $seed_hmm_dir does not exists create it.
#print "Creating $seed_hmm_dir\n" unless -e $seed_hmm_dir;
#`mkdir -p $seed_hmm_dir`         unless -e $seed_hmm_dir;
##if the $aln_stock_dir does not exists create it.
#print "Creating $aln_stock_dir\n" unless -e $aln_stock_dir;
#`mkdir -p $aln_stock_dir`         unless -e $aln_stock_dir;

if ( $type eq 'new' ) {
	for ( my $i = $low_bound; $i <= $high_bound; $i = $i + $step ) {
		print "processing family : $i\t$directory/$i$suffix\n";
		my $fc = get_family_construction_id($i);
		if ( -e "$directory/$i$suffix" ) {
			if ( -e "$rep_dir/$i.rep$suffix" ) {
				##copy the representatives to the repository
				`cp $rep_dir/$i.rep$suffix $output_dir/FC_$fc/seqs_reps/$i.faa`;
				##Align the representatives
				`muscle -in $rep_dir/$i.rep$suffix -out $output_dir/FC_$fc/alns_seed/$i.aln` unless -e "$output_dir/FC_$fc/alns_seed/$i.aln";
				##Build HMM of representatives Set this as seed HMM.
				`hmmbuild --informat afa $output_dir/FC_$fc/hmms_seed/$i.hmm $output_dir/FC_$fc/alns_seed/$i.aln` unless -e "$output_dir/FC_$fc/hmms_seed/$i.hmm";
				##Align non representatives to HMM
				`hmmalign --outformat afa --trim --amino -o $output_dir/FC_$fc/alns_full/$i.aln $output_dir/FC_$fc/hmms_seed/$i.hmm $directory/$i$suffix`;
			} else {
				##align the fasta file
				`muscle -in $directory/$i$suffix -out $output_dir/FC_$fc/alns_full/$i.aln` unless -e "$output_dir/FC_$fc/alns_full/$i.aln";
				##copy the alignment to the seed directory
				`cp $output_dir/FC_$fc/alns_full/$i.aln $output_dir/FC_$fc/alns_seed/$i.aln`;
			}
			##transform the fasta aln to a stockholm aln
			fasta_aln_to_stockholm( file_in => "$output_dir/FC_$fc/alns_full/$i.aln", file_out => "$output_dir/FC_$fc/alns_stock_full/$i.stock" );
			fasta_aln_to_stockholm( file_in => "$output_dir/FC_$fc/alns_seed/$i.aln", file_out => "$output_dir/FC_$fc/alns_stock_seed/$i.stock" );
			##Build full HMM
			`hmmbuild --informat afa $output_dir/FC_$fc/hmms_full/$i.hmm $output_dir/FC_$fc/alns_full/$i.aln` unless -e "$output_dir/FC_$fc/hmms_full/$i.hmm";
			`cp $output_dir/FC_$fc/hmms_full/$i.hmm $output_dir/FC_$fc/hmms_seed/$i.hmm` unless -e "$output_dir/FC_$fc/hmms_seed/$i.hmm" ;
			if ( -e "$output_dir/FC_$fc/hmms_full/$i.hmm" && !-e "$output_dir/FC_$fc/hmms_full/$i.hmm.gz" ) {
				##zip the file
				`gzip $output_dir/FC_$fc/hmms_full/$i.hmm`;
			}
			if ( -e "$output_dir/FC_$fc/hmms_seed/$i.hmm" && !-e "$output_dir/FC_$fc/hmms_seed/$i.hmm.gz" ) {
				##zip the file
				`gzip $output_dir/FC_$fc/hmms_seed/$i.hmm`;
			}
			##Tree of the alignment
			`FastTree $output_dir/FC_$fc/alns_full/$i.aln > $output_dir/FC_$fc/trees_full/$i.tree` unless -e "$output_dir/FC_$fc/trees_full/$i.tree";
		}
	}
} elsif ( $type eq 'old' ) {
	for ( my $i = $low_bound; $i <= $high_bound; $i = $i + $step ) {
		my $fc = get_family_construction_id($i);
		##unzip the hmm if it is zipped
		`gunzip $ref_hmm_dir/FC_$fc/hmms_full/$i.hmm.gz` if -e "$ref_hmm_dir/FC_$fc/hmms_full/$i.hmm.gz";
		##align sequences  Might us --trim later
		`hmmalign --outformat afa --amino -o $output_dir/FC_$fc/alns_full/$i.aln --mapali $ref_hmm_dir/FC_$fc/alns_stock_full/$i.stock $ref_hmm_dir/FC_$fc/hmms_full/$i.hmm $directory/$i$suffix`
		  if -e "$directory/$i$suffix";
		##zip the hmm after it has been used
		`gzip $ref_hmm_dir/FC_$fc/hmms_full/$i.hmm` if -e "$ref_hmm_dir/FC_$fc/hmms_full/$i.hmm";
		## transform the fasta aln into a stockholm aln
		fasta_aln_to_stockholm( file_in => "$output_dir/FC_$fc/alns_full/$i.aln", file_out => "$output_dir/FC_$fc/alns_stock_full/$i.stock" ) unless -s "$output_dir/FC_$fc/alns_stock_full/$i.stock";
		##HMM of the alignment file
		`hmmbuild --informat afa $output_dir/FC_$fc/hmms_full/$i.hmm $output_dir/FC_$fc/alns_full/$i.aln` unless -s "$output_dir/FC_$fc/hmms_full/$i.hmm";
		##compress the HMM file
		`gzip $output_dir/FC_$fc/hmms_full/$i.hmm` if -e "$output_dir/FC_$fc/hmms_full/$i.hmm";
		##Tree of the alignment
		`FastTree $output_dir/FC_$fc/alns_full/$i.aln > $output_dir/FC_$fc/trees_full/$i.tree` unless -s "$output_dir/FC_$fc/trees_full/$i.tree";
	}
} else {
	die "Command $type not recognized\n";
}


sub get_family_construction_id{
	my $famid = shift;
	foreach my $fcid (keys %fc_range){
		return $fcid if $fc_range{$fcid}[0] <= $famid && $fc_range{$fcid}[1] >= $famid;
	}
}

sub get_family_size_by_sql {
	my %args   = @_;
	my $fam_id = $args{famid};

	#my $family = $analysis->MRC::DB::get_family_from_famid($fam_id);
	#return $family>get_column('size');
}

sub get_family_size_by_grep {
	my %args = @_;
	my $fam  = $args{family};
	return `grep -c '>' $fam`;
}

sub fasta_aln_to_stockholm {
	my %args     = @_;
	my $file_in  = $args{file_in};
	my $file_out = $args{file_out};
	my $io_obj   = Bio::AlignIO->new( -file => $file_in, -format => 'fasta' );
	my $out_obj  = Bio::AlignIO->new( -file => ">$file_out", -format => 'stockholm' );
	while ( my $aln = $io_obj->next_aln() ) {
		$out_obj->write_aln($aln);
	}
	return 1;
}
