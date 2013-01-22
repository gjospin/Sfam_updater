#1/usr/bin/perl -w
package Sfam_updater;
use strict;
use Sfam_updater::launch_sifting;

#use Sfam_updater::DB_op;
#use Sfam_updater::download_jgi_data;
use Getopt::Long;

#gather options for what to skip. If nothing is specified, run everything.
my (
	 $username,              $password,            $skip_download,            $skip_index,           $skip_sifting,
	 $skip_DB_inserts,       $total_files_to_sift, $core_file_to_sift,        $total_files_to_index, $core_file_to_index,
	 $genome_download_dir,   $genome_file,         $tmp_data,                 $sifting_lastal_db,    $lastal_results_core,
	 $lastal_leftovers_core, $CDS_length_hash_ref, $lastal_new_familymembers, $hmm_file_core,        $number_of_hmms,
	 $hmmsearch_results_core
);

#my $DB_pointer  = "DBI:mysql:Sfams:lighthouse.ucsf.edu";
my $DB_pointer = "DBI:mysql:Sfams";
GetOptions(
	"skip-download"   => \$skip_download,
	"skip-index"      => \$skip_index,
	"skip-sifting"    => \$skip_sifting,
	"skip-DB-inserts" => \$skip_DB_inserts,
	"u=s"             => \$username,
	"p=s"             => \$password,
	"download-dir=s"  => \$genome_download_dir,
	"genome-file=s"   => \$genome_file,
	"data-dir=s"      => \$tmp_data,

);
$password = <>;
chomp($password);

# Read IMG excel spreadsheet and download new genomes
# Insert new genomes and genes information in the DB
# Extract (and index) new CDS

unless ($skip_download) {

	#	my $new_genome_oids_ref = Sfam_updater::DB_op::add_new_genomes(
	#																	username    => $username,
	#																	password    => $password,
	#																	db          => $DB_pointer,
	#																	genome_file => $genome_file
	#	);
	#	$new_genome_oids_ref = Sfam_updater::download_jgi_data::download_data(
	#																		   username         => $username,
	#																		   password         => $password,
	#																		   db               => $DB_pointer,
	#																		   db_master        => $genome_download_dir,
	#																		   genome_oid_array => $new_genome_oids_ref,
	#	);
	#	Sfam_updater::DB_op::add_new_genes(
	#										password         => $password,
	#										db               => $DB_pointer,
	#										file_repo        => $genome_download_dir,
	#										genome_oid_array => $new_genome_oids_ref,
	#										new_cds_dump_dir => $tmp_data,
	#										db_master        => $genome_download_dir,
	#	);

}

unless ($skip_index) {

	#( $total_files_to_sift, $core_file_to_sift ) = Sfam_updater::DB_op::gather_CDS(
	#																				output_dir => "/home/gjospin/proteinFamilies/Sfam_updater/testdir",
	#																				old        => 0,
	#																				db         => $DB_pointer,
	#																				username   => $username,
	#																				password   => $password,
	#																				fragmented => 1
	#);
	#	( $total_files_to_index, $core_file_to_index ) = Sfam_updater::DB_op::gather_CDS(
	#																					  output_dir => $tmp_data,
	#																					  old        => 1,
	#																					  db         => $DB_pointer,
	#																					  username   => $username,
	#																					  password   => $password,
	#																					  fragmented => 0
	#	);

	# Sift new CDS using Lastal
	$sifting_lastal_db = Sfam_updater::launch_sifting::index_familymembers(    file_to_index => "$core_file_to_index.fasta",
																			output_dir    => $tmp_data."/ref_lastal_db" );
}
$sifting_lastal_db   = "/share/eisen-z2/gjospin/Sfam_updater/test_dir/ref_lastal_db/sfams";
$core_file_to_index = "/share/eisen-z2/gjospin/Sfam_updater/test_dir/familymembers";
$core_file_to_sift   = "/share/eisen-z2/gjospin/Sfam_updater/test_dir/newCDS";
my $data_repo = "/var/opt/iseem/proteinFamilies/HMMs";
unless ($skip_sifting) {

	$lastal_results_core = Sfam_updater::launch_sifting::launch_lastal(
												db            => $sifting_lastal_db,
												filename_core => $core_file_to_sift,
												arguments     => "-e100 -m100 -l3 -f0",
												output_dir    => "$tmp_data/lastal_sift_output",
												error_dir     => "$tmp_data/lastal_sift_err"
	);
	$lastal_results_core = "/share/eisen-z2/gjospin/Sfam_updater/test_dir/lastal_sift_output/last_sift.ouput";

	#$lastal_results_core = "/share/eisen-z2/gjospin/Sfam_updater/test_dir/lastal_sift_output/last_sift.ouput";
	$CDS_length_hash_ref = Sfam_updater::launch_sifting::compile_sequence_length_hash( file1 => "$core_file_to_index", file2 => " $core_file_to_sift" );
	$lastal_new_familymembers = Sfam_updater::launch_sifting::parse_lastal(
		filename_core => $lastal_results_core,
		output_dir    => "$tmp_data/lastal_sift_output",
		length_hash   => $CDS_length_hash_ref

	);

	#Sfam_updater::DB_op::add_new_family_members(
	#											 family_members_file => $lastal_new_familymembers,
	#											 db                  => $DB_pointer,
	#											 username            => $username,
	#											 password            => $password,
	#);
	$lastal_leftovers_core = Sfam_updater::launch_sifting::extract_lastal_leftovers(
																					 family_members_files => $lastal_new_familymembers,
																					 new_cds_core         => $core_file_to_sift,
																					 output_dir           => $tmp_data."/lastal"
	);

	# Break HMMs into manageable size files
	( $number_of_hmms, $hmm_file_core ) = Sfam_updater::launch_sifting::index_hmms(
		db         => $DB_pointer,
		username   => $username,
		password   => $password,
		output_dir => $tmp_data."/hmmsearch",
		file_size  => 2000,
		repo       => $data_repo

	);

	# Launch HMMsearch of each lastal_leftover file again all HMM files
	$hmmsearch_results_core = Sfam_updater::launch_sifting::launch_hmmsearch(
																			  hmmfiles   => $hmm_file_core,
																			  seq_files  => $lastal_leftovers_core,
																			  arguments  => "-E 1e-5 -Z $number_of_hmms --domtblout ",
																			  output_dir => "$tmp_data/hmmsearch_sift_output",
																			  error_dir  => "$tmp_data/hmmsearch_sift_err"
	);
}

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
