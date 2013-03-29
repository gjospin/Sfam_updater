#!/usr/bin/perl -w
package Sfam_updater;
use strict;
use Sfam_updater::launch_sifting;

use Sfam_updater::DB_op;
use Sfam_updater::download_jgi_data;
use Getopt::Long;

#dependencies in Sfam_updater include the MRC perl modules (github.com/sharpton/MRC)
#install them as follows:
#% git clone git@github.com:sharpton/MRC.git
#and add to your PERL5LIB path

#MUST HAVE MCL INSTALLED ON THE LOCAL MACHINE!
#MUST HAVE HMMSEARCH, BLAST, AND FASTTREE INSTALLED ON THE REMOTE MACHINE!

my ( $tmp_data );
my $threads            = 2; #default hmmsearch_thread
my $rep_threshold      = 250;
my $path_to_sfams_repo = "./";
GetOptions(
    "data-dir|d=s"      => \$tmp_data,
    "threads|t=i"       => \$threads,
    "rep-thresh|r=i"    => \$rep_threshold,
);

print "Processing data in $tmp_data\n";

#REP PICKING NEEDS TO BE DONE ON REMOTE SERVER
my $representatives_dir = Sfam_updater::DB_op::prep_families_for_representative_picking_nosql_nopremcl(
    output_directory => $tmp_data,
    rep_threshold    => $rep_threshold,	
    blast_args       => "-outfmt 6",
    threads          => $threads,
    fasta_stem       => ".fa",
    );

#my $representatives_dir = $tmp_data."/new_fams/representatives";

die;

Sfam_updater::launch_sifting::fine_tune_representative_picking( 
    reps_dir        => $representatives_dir,
    get_link_path   => $path_to_sfams_repo . "/get_link_by_list.pl",
    mcl_redunt_path => $path_to_sfams_repo . "/mcl_redunt_reduce.pl",
    rep_threshold   => $rep_threshold,
    );

Sfam_updater::DB_op::generate_representative_fasta_nosql(
    representative_dir => $representatives_dir,
    output_dir         => $tmp_data,
    fasta_stem         => ".fa",
    );

#De novo alignment for $tmp_data/new_fams files
my $new_fam_dir = Sfam_updater::launch_sifting::build_aln_hmm_trees(
    directory  => $tmp_data,
    repo       => $tmp_data,
    total_jobs => 200,
    type       => 'new',
    output     => $tmp_data."new_fams",
    error      => $tmp_data."aln_hmm_trees_new",
    );
