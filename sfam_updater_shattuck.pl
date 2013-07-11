#!/usr/bin/perl -w
package Sfam_updater;
use strict;
use Sfam_updater::launch_sifting;

use Sfam_updater::DB_op;
use Sfam_updater::download_jgi_data;
use Getopt::Long;

#MUST HAVE MCL INSTALLED ON THE LOCAL MACHINE!
#MUST HAVE HMMSEARCH, BLAST, FASTTREE, AND MCL INSTALLED ON THE REMOTE MACHINE!

#gather options for what to skip. If nothing is specified, run everything.
my (
    $username,               $password,             $skip_download,            $skip_index,           $skip_sifting,
    $skip_DB_inserts,        $total_files_to_sift,  $core_file_to_sift,        $total_files_to_index, $core_file_to_index,
    $genome_download_dir,    $genome_file,          $tmp_data,                 $sifting_lastal_db,    $lastal_results_core,
    $lastal_leftovers_core,  $CDS_length_hash_ref,  $lastal_new_familymembers, $hmm_file_core,        $number_of_hmms,
    $hmmsearch_results_core, $total_files_to_blast, $core_file_to_blast,       $author,               $description,
    $name,                   $data_repo,            $rep_threshold,            $skip_hmm_index,       
);
my $threads            = 2;      #default hmmsearch_thread
my $path_to_sfams_repo = "./";

my $family_construction_id;

#my $DB_pointer  = "DBI:mysql:Sfams:lighthouse.ucsf.edu";#my $DB_pointer = "DBI:mysql:SFams";
my $DB_pointer = "DBI:mysql:SFams_dev";
#my $DB_pointer = "DBI:mysql:SFams_MetaHIT";

GetOptions(
    "skip-download"   => \$skip_download,
    "skip-index"      => \$skip_index,
    "skip-sifting"    => \$skip_sifting,
    "skip-DB-inserts" => \$skip_DB_inserts,
    "skip-hmm-index"  => \$skip_hmm_index,
    "u=s"             => \$username,
    "p=s"             => \$password,
    "db=s"            => \$DB_pointer,
    "download-dir=s"  => \$genome_download_dir,
    "genome-file=s"   => \$genome_file,
    "data-dir=s"      => \$tmp_data,
    "ref-db=s"        => \$data_repo,             #location of the ref hmms
    "threads=i"       => \$threads,
    "author=s"        => \$author,
    "description=s"   => \$description,
    "name=s"          => \$name,
    "rep-thresh=i"    => \$rep_threshold,
    "fci:i"           => \$family_construction_id, #if not defined, we autoset a new one.
    );

die "Please provide a name for the new family construction using --name\n"                   unless defined($name);
die "Please provide an Author name for the new family construction using --author\n"         unless defined($author);
die "Please provide a description for the new family construction using --description\n"     unless defined($description);
die "Please provide a username for the MySQL database (Needs insert priviledges) using -u\n" unless defined($username);
die "Please provide a Database pointer for the MySQL database to use using --db\n"           unless defined($DB_pointer);
die "Please provide the location of the reference HMM database repository using --ref-db\n"  unless defined($data_repo) || ($skip_index);
die "Please define a representative selection threshold by using --rep-thresh\n"             unless defined($rep_threshold);

if ( !defined($password) ) {
	print "Enter MySQL password :\n";
	$password = <>;
	chomp($password);
}

#first thing, create a familyconstruction row
if( !defined( $family_construction_id ) ){
    $family_construction_id = Sfam_updater::DB_op::insert_fc(
	username    => $username,
	password    => $password,
	db          => $DB_pointer,
	author      => $author,
	description => $description,
	name        => $name
	);
}

# Read IMG excel spreadsheet and download new genomes
# Insert new genomes and genes information in the DB
# Extract (and index) new CDS
prepare_data_directory( fc_id => $family_construction_id, output => $tmp_data );
unless ($skip_download) {
    die "Please provide a directory path for the flat file repository that we'll dump downloaded data to using --download-dir\n"
	unless defined($genome_download_dir);
    my $new_genome_oids_ref = Sfam_updater::DB_op::add_new_genomes(
	username    => $username,
	password    => $password,
	db          => $DB_pointer,
	genome_file => $genome_file
	);
    $new_genome_oids_ref = Sfam_updater::download_jgi_data::download_data(
	username         => $username,
	password         => $password,
	db               => $DB_pointer,
	db_master        => $genome_download_dir,
	genome_oid_array => $new_genome_oids_ref,
	);
    Sfam_updater::DB_op::add_new_genes(
	password         => $password,
	db               => $DB_pointer,
	file_repo        => $genome_download_dir,
	genome_oid_array => $new_genome_oids_ref,
	new_cds_dump_dir => $tmp_data,
	db_master        => $genome_download_dir,
	);
}

unless ($skip_index) {
    $core_file_to_sift = Sfam_updater::DB_op::gather_CDS(
	output_dir => $tmp_data,
	old        => 0,
	db         => $DB_pointer,
	username   => $username,
	password   => $password,
	fragmented => 1
	);
    
    #$core_file_to_index = Sfam_updater::DB_op::gather_CDS(
    #																					  output_dir => $tmp_data,
    #																					  old        => 1,
    #																					  db         => $DB_pointer,
    #																					  username   => $username,
    #																					  password   => $password,
    #																					  fragmented => 0
    #	);
    
    # Sift new CDS using Lastal
    #$sifting_lastal_db = Sfam_updater::launch_sifting::index_familymembers(    file_to_index => "$core_file_to_index.fasta",
    #																	output_dir    => $tmp_data."/ref_lastal_db" );
}


my $new_family_members;
unless ($skip_sifting) {
    unless( $skip_hmm_index ){
	( $number_of_hmms, $hmm_file_core ) = Sfam_updater::launch_sifting::index_hmms(
	    db         => $DB_pointer,
	    username   => $username,
	    password   => $password,
	    output_dir => $tmp_data."/hmmsearch",
	    file_size  => 2000,
	    repo       => $data_repo
	    );
    }
    
    #$lastal_leftovers_core = "/share/eisen-z2/gjospin/Sfam_updater/test_dir/lastal_leftovers.fasta";
    
    # Launch HMMsearch of each lastal_leftover file again all HMM files
    #assume that the cluster can read the files directly from our db server
    #for shattuck, need to have the data and script in the folowing location on shattuck
    #/mnt/data/work/pollardlab
    #chef will see this as the following location
    #/pollard/shattuck0
    #$core_file_to_sift = "/mnt/data/work/pollardlab/sharpton/Sfam_updater/test/sfams/newCDS_1.fasta";
    #$hmm_file_core     = "/mnt/data/work/pollardlab/sharpton/Sfam_updater/test/sfams/hmmsearch/HMMs";
    #$number_of_hmms    = 10;
    
    #need to go into the script creation process and fix the run-time settings, which are currently set to short.q for testing.
    if( $skip_index ){
	$core_file_to_sift = "${tmp_data}/newCDS_1.fasta"; 
    }
    if( $skip_hmm_index ){
	$hmm_file_core    = "${tmp_data}/hmmsearch/HMMs";
#	$number_of_hmms   = count_hmms( $hmm_file_core );
	$number_of_hmms   = 436351;
    }
    unless( $skip_sifting ){
	$hmmsearch_results_core = Sfam_updater::launch_sifting::launch_hmmsearch(
	    hmmfiles   => $hmm_file_core,
	    seq_files  => $core_file_to_sift,
	    arguments  => "-E 1e-5 -Z $number_of_hmms ",
	    output_dir => "$tmp_data/hmmsearch_sift_output",
	    error_dir  => "$tmp_data/hmmsearch_sift_err",
	    threads    => $threads,
	    machine    => "chef",
	    pull_only  => 0,  #set this if you only want to grab remote results, in case that handler script died during cluster's run
	    );
	
	#    $hmmsearch_results_core = "/mnt/data/work/pollardlab/sharpton/Sfam_updater/test/sfams/hmmsearch_sift_output/hmmsearch_sift.ouput";
	$new_family_members = Sfam_updater::launch_sifting::parse_hmmsearch(       
	    filename_core => $hmmsearch_results_core,
	    output_dir    => "$tmp_data/hmmsearch_sift_output", 
	    );
    }
}
if( $skip_sifting ){
    $new_family_members =  $tmp_data . "/hmmsearch_sift_output/hmmsearch_newCDS_to_fam.map";
}
Sfam_updater::DB_op::insert_familymembers(
    input    => $new_family_members, #table that maps gene_oid to famid, tab delimited
    db       => $DB_pointer,
    username => $username,
    password => $password,
    output   => $tmp_data."/old_fams",
    );
die;
$core_file_to_blast = Sfam_updater::DB_op::gather_CDS(
    output_dir => "$tmp_data/blast_input",
    old        => 0,
    db         => $DB_pointer,
    username   => $username,
    password   => $password,
    fragmented => 1,
    spl_size   => 200,
    );

my $total_seqs_to_blast = Sfam_updater::DB_op::count_all_CDS(
    old      => 0,
    db       => $DB_pointer,
    username => $username,
    password => $password,
    );

#my $total_seqs_to_blast = 279555;

die; #for testing...

#need to go into the script creation process and fix the run-time settings, which are currently set to short.q for testing.
#need to fix data transfer
my $blast_results_core = Sfam_updater::launch_sifting::launch_blast(
    output_dir             => "$tmp_data/blast_results",
    blast_input_files_core => $core_file_to_blast,
    arguments => "-outfmt 6 -searchsp $total_seqs_to_blast -num_alignments 100000 ",
    error_dir => "$tmp_data/blast_err",
    threads   => $threads,
    machine   => "chef",
    );

#    my $blast_results_core = "/share/eisen-z2/gjospin/Sfam_updater/test_dir/blast_results/blast_output";
#my $blast_results_core = "/mnt/data/work/pollardlab/sharpton/Sfam_updater/test/sfams/blast_results/blast_output";
#$core_file_to_blast = "/mnt/data/work/pollardlab/sharpton/Sfam_updater/test/sfams/blast_input/newCDS";

my $blast_seqs_lengths = Sfam_updater::launch_sifting::compile_sequence_length_hash( file2 => $core_file_to_blast );

my $mcl_file = Sfam_updater::launch_sifting::parse_blast(
    output_dir         => "$tmp_data/MCL",
    blast_results_core => $blast_results_core,
    coverage           => 0.8,
    evalue             => "1e-10",
    seq_lengths        => $blast_seqs_lengths,
    );

#might need to think more about the memory and run time reqirements associated with this step
my $mcl_output_file = Sfam_updater::launch_sifting::launch_mcl(
    output_dir => "$tmp_data/MCL",
    
    #	    queue      => "-l jumbo -l h_vmem=50G",
    queue      => "-l mem_free=1G",
    input      => $mcl_file,
    mcl_params => "-I 2.0",
    error_dir  => "$tmp_data/MCL",
    threads    => 2,
    machine    => "chef",
    );

#my $mcl_output_file             = "/home/gjospin/proteinFamilies/Sfam_updater/merlot_test/test_dir/MCL/mcl_output.mcl";
#my $new_family_IDs_mapping_file = "/home/gjospin/proteinFamilies/Sfam_updater/merlot_test/test_dir/MCL/mcl_newCDS_to_fam.map";
#my $mcl_output_file = "/mnt/data/work/pollardlab/sharpton/Sfam_updater/test/sfams/MCL/mcl_output.mcl";
my $new_family_IDs_mapping_file = Sfam_updater::launch_sifting::parse_mcl(
    db              => $DB_pointer,
    username        => $username,
    password        => $password,
    mcl_output_file => $mcl_output_file,
    output_dir      => $tmp_data."/MCL",
    );
Sfam_updater::DB_op::print_fc_range(
    db         => $DB_pointer,
    username   => $username,
    password   => $password,
    output_dir => $tmp_data,
    ); 
Sfam_updater::DB_op::insert_familymembers(
    input    => $new_family_IDs_mapping_file,
    db       => $DB_pointer,
    username => $username,
    password => $password,
    output   => $tmp_data."/new_fams",
    );

# Seqs have been split by fam_ids.  Need to align
#De novo alignment for $tmp_data/new_fams files

#$family_construction_id = 65;
#my $mcl_file = "/home/gjospin/proteinFamilies/Sfam_updater/merlot_test/test_dir/MCL_input/mcl_input.abc";
#my $mcl_file = "/mnt/data/work/pollardlab/sharpton/Sfam_updater/test/sfams/MCL/mcl_input.abc";

#REP PICKING NEEDS TO BE DONE ON REMOTE SERVER
my $representatives_dir = Sfam_updater::DB_op::prep_families_for_representative_picking(
    db               => $DB_pointer,
    username         => $username,
    password         => $password,
    output_directory => $tmp_data."/new_fams",
    fc_id            => $family_construction_id,
    mcl_input        => $mcl_file,
    rep_threshold    => $rep_threshold,
    );

#my $representatives_dir = $tmp_data."/new_fams/representatives";

Sfam_updater::launch_sifting::fine_tune_representative_picking(
    reps_dir        => $representatives_dir,
    get_link_path   => $path_to_sfams_repo."/get_link_by_list.pl",
    mcl_redunt_path => $path_to_sfams_repo."/mcl_redunt_reduce.pl",
    rep_threshold   => $rep_threshold,
    );
Sfam_updater::DB_op::generate_representative_fasta(
    representative_dir => $representatives_dir,
    output_dir         => $representatives_dir,
    db                 => $DB_pointer,
    username           => $username,
    password           => $password,
    );

#this isn't working yet...
my $old_fam_dir = Sfam_updater::launch_sifting::build_aln_hmm_trees(
    directory  => $tmp_data."/old_fams",
    repo       => $data_repo,
    total_jobs => 200,
    type       => 'old',
    output     => $tmp_data."old_fams",
    error      => $tmp_data."aln_hmm_trees_old",
    machine    => "chef",
    );

#De novo alignment for $tmp_data/new_fams files
my $new_fam_dir = Sfam_updater::launch_sifting::build_aln_hmm_trees(
    directory  => $tmp_data."/new_fams",
    repo       => $data_repo,
    total_jobs => 200,
    type       => 'new',
    output     => $tmp_data."new_fams",
    error      => $tmp_data."aln_hmm_trees_new",
    );

## Insert tree into DB
#	my $new_fam_dir = "/home/gjospin/proteinFamilies/Sfam_updater/merlot_dir/test_dir/new_fams";
#my $new_fam_dir = "/home/gjospin/proteinFamilies/Sfam_updater/merlot_test/test_dir/new_fams";
#>>>>>>> ee06a0ca8ed336755e0bdd20fb85f88979d38c9e
Sfam_updater::DB_op::insert_trees_hmms_alignments(
    directory           => $new_fam_dir,
    db                  => $DB_pointer,
    username            => $username,
    password            => $password,
    tree_desc           => "Tree build for new families for fci $family_construction_id",
    tree_type           => "alltree",
    tree_path           => "/trees",
    aln_path            => "/alignments",
    seed_alignment_path => "/seed_alignments",
    hmm_path            => "/HMMs",
    );

#Sfam_updater::DB_op::insert_trees(
#								   directory => $old_fam_dir,
#								   db        => $DB_pointer,
#								   username  => $username,
#								   password  => $password,
#								   tree_desc => "Tree build for old families for fci $family_construction_id",
#								   tree_type => "alltree",
#								   tree_path => "/trees";
    #);
    
    # Insert alignment into DB
    #	Sfam_updater::DB_op::insert_alignments(
    #	    directory => $new_fam_dir,
    #	    db        => $DB_pointer,
    #	    username  => $username,
    #	    password  => $password,
    #	    type=> "new",
    #	);
    #      Sfam_updater::DB_op::insert_alignments(
    #            directory => $old_fam_dir,
    #            db        => $DB_pointer,
    #            username  => $username,
    #            password  => $password,
    #	  type => "old",
    #	  );
    
    # Insert hmm into DB
    # Insert seed hmm into DB for new families
    
    # package update / Release.
    ## move everything from the old DB into a subdirectory.
    ## populate "current" DB release directory with new data.

exit;

# Parse lastal results
# Add new familymembers information
# Sift left over new sequences using HMMsearch
# Parse HMMsearch results
# Add new familymembers information
# Extract ALL CDS that are not familymembers + lastal Index
# Launch All versus All using Lastal
# Parse lastal results + prep MCL input
# Run MCL on a large memory node (If needed)
# Insert new families into DB
# Add new familymembers information
# Launch build of new alignments
# Launch build of new HMMs
# Launch build of new Trees
# Launch build of updated alignments
# Launch build of updated HMMs
# Launch build of updated Trees

sub prepare_data_directory {
	my %args      = @_;
	my $max_famid = $args{fc_id};
	my $output   = $args{output};
	for ( my $i = $max_famid; $i >= 0; $i-- ) {
		my $new_repo = $output."/FC_$i";
		my $aln = $new_repo."/alns_full";
		`mkdir -p $aln` unless -e $aln;
		my $aln_stock = $new_repo."/alns_stock_full";
		`mkdir -p $aln_stock` unless -e $aln_stock;
		my $aln_stock_seed = $new_repo."/alns_stock_seed";
		`mkdir -p $aln_stock_seed` unless -e $aln_stock_seed;
		my $seed_aln = $new_repo."/alns_seed";
		`mkdir -p $seed_aln` unless -e $seed_aln;
		my $code = $new_repo."/code";
		`mkdir -p $code`unless -e $code;
		my $function = $new_repo."/function";
		`mkdir -p $function`unless -e $function;
		my $hmms_full = $new_repo."/hmms_full";
		`mkdir -p $hmms_full`unless -e $hmms_full;
		my $hmms_seed = $new_repo."/hmms_seed";
		`mkdir -p $hmms_seed`unless -e $hmms_seed;
		my $mysql_db_dump = $new_repo."/mysql_db_dump";
		`mkdir -p $mysql_db_dump`unless -e $mysql_db_dump;
		my $seqs_all = $new_repo."/seqs_all";
		`mkdir -p $seqs_all`unless -e $seqs_all;
		my $seqs_reps = $new_repo."/seqs_reps";
		`mkdir -p $seqs_reps`unless -e $seqs_reps;
		my $stats = $new_repo."/stats";
		`mkdir -p $stats`unless -e $stats;
		my $trees_full = $new_repo."/trees_full";
		`mkdir -p $trees_full`unless -e $trees_full;
		my $precision = $stats."/precision_recall";
		`mkdir -p $precision`unless -e $precision;
		my $clans = $stats."/clans";
		`mkdir -p $clans`unless -e $clans;
		my $kegg = $function."/kegg";
		`mkdir -p $kegg`unless -e $kegg;
		my $interpro = $function."/interpro";
		`mkdir -p $interpro`unless -e $interpro;
	}
}

sub count_hmms{
    my $hmm_file_core  = shift;
    my $number_of_hmms = 0;
    my @files = glob( "${hmm_file_core}*" );
    foreach my $file( @files ){
	my $count = `grep -c "NAME" $file`;
	chomp $count;
	$number_of_hmms += $count;
    }
    print "Found $number_of_hmms HMMs\n";
    return $number_of_hmms;
}
