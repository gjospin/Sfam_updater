#!/usr/bin/perl -w
package Sfam_updater::download_jgi_data;
use strict;
use Getopt::Long;
use MRC;
use MRC::DB;
use DBI;

# Code written by Thomas Sharpton
# Code integrated in the Sfam_updater by Guillaume Jospin

sub download_data {
	my %args            = @_;
	my $username        = $args{username};
	my $password        = $args{password};
	my $lookup_seq_file = $args{lookup_seq_file};
	my $ffdb_master     = $args{db_master};
	my $db_pointer      = $args{db};
	my $analysis        = MRC->new();
	my $genome_oid_ref  = $args{genome_oid_array};
	my %genome_oids     = %{$genome_oid_ref};

	#setting up the DB connection
	$analysis->set_dbi_connection($db_pointer);
	$analysis->set_username($username);
	$analysis->set_password($password);
	$analysis->build_schema();

	#print "Before connection\n$db_pointer\t$username\t$password\n";;

	#my $DB = DBI->connect( "$db_pointer", "$username", "$password" ) or die "Couldn't connect to database : ".DBI->errstr;
	#print "After connection\n";
	foreach my $taxon_oid ( keys %genome_oids ) {

		#my $statement = $DB->prepare('SELECT directory FROM genomes WHERE genome_oid = ?');
		#$statement->execute($taxon_oid);
		#my $taxon_dir = $statement->fetchrow();
		my $genome = $analysis->MRC::DB::get_genome_from_taxon_oid($taxon_oid);
		my $taxon_dir;
		$taxon_dir = $genome->directory();
		my $taxon_full_dir = $ffdb_master.$taxon_dir."/";
		system("mkdir -p $taxon_full_dir");
		print "Downloading $taxon_dir\n";
		get_proteins( $taxon_oid, $taxon_full_dir );
		get_genome( $taxon_oid, $taxon_full_dir );
		get_genes( $taxon_oid, $taxon_full_dir );
		get_gene_info( $taxon_oid, $taxon_full_dir );
		get_intergenic_seqs( $taxon_oid, $taxon_full_dir );
	}
	return \%genome_oids;
}

sub get_proteins {
	my ( $taxon_oid, $taxon_full_dir ) = @_;
	my $output = $taxon_full_dir.$taxon_oid.".faa";

	#return if -e $output.".gz";
	my $url = "\"http://img.jgi.doe.gov/cgi-bin/w/main.cgi?section=TaxonDetail&downloadTaxonFaaFile=1&_noHeader=1&taxon_oid=".$taxon_oid."\"";
	if ( -e $output.".gz" ) {
		print STDERR "$output.gz exists \n";
		if ( -e $output ) {
			print STDERR "Old output exists removing\n";
			`rm $output`;
		}
		print STDERR "Unzipping $output.gz\n";
		system("gunzip $output.gz");
		print STDERR "Done Unzipping $output.gz\n";
		if ( -z $output ) {
			print STDERR "File is empty, downloading again\n";
			download_url( $url, $output );
		}
	} else {
		download_url( $url, $output );
	}
	print STDERR "Zipping $output\n";
	`gzip $output`;
	print STDERR "Done Zipping $output\n";
	if ( !( -e $output.".gz" ) ) {
		warn("Couldn't download $output from:\n$url\n");
	}
}

sub get_genome {
	my ( $taxon_oid, $taxon_full_dir ) = @_;
	my $output = $taxon_full_dir.$taxon_oid.".ffa";

	#return if -e $output."gz";
	my $url = "\"http://img.jgi.doe.gov/cgi-bin/w/main.cgi?section=TaxonDetail&downloadTaxonFnaFile=1&_noHeader=1&taxon_oid=".$taxon_oid."\"";
	if ( -e $output.".gz" ) {
		if ( -e $output ) {
			`rm $output`;
		}
		system("gunzip $output.gz");
		if ( -z $output ) {
			download_url( $url, $output );
		}
	} else {
		download_url( $url, $output );
	}
	`gzip $output`;
	if ( !( -e $output.".gz" ) ) {
		warn("Couldn't download $output from:\n$url\n");
	}
}

sub get_genes {
	my ( $taxon_oid, $taxon_full_dir ) = @_;
	my $output = $taxon_full_dir.$taxon_oid.".ffn";

	#return if -e $output."gz";
	my $url = "\"http://img.jgi.doe.gov/cgi-bin/w/main.cgi?section=TaxonDetail&downloadTaxonGenesFnaFile=1&_noHeader=1&taxon_oid=".$taxon_oid."\"";
	if ( -e $output.".gz" ) {
		if ( -e $output ) {
			`rm $output`;
		}
		system("gunzip $output.gz");
		if ( -z $output ) {
			download_url( $url, $output );
		}
	} else {
		download_url( $url, $output );
	}
	`gzip $output`;
	if ( !( -e $output.".gz" ) ) {
		warn("Couldn't download $output from:\n$url\n");
	}
}

sub get_gene_info {
	my ( $taxon_oid, $taxon_full_dir ) = @_;
	my $output = $taxon_full_dir.$taxon_oid.".genes_info";

	#return if -e $output."gz";
	my $url = "\"http://img.jgi.doe.gov/cgi-bin/w/main.cgi?section=TaxonDetail&downloadTaxonGenesFile=1&_noHeader=1&taxon_oid=".$taxon_oid."\"";
	if ( -e $output.".gz" ) {
		if ( -e $output ) {
			print STDERR "Removing output\n";
			`rm $output`;
		}
		system("gunzip $output.gz");
		if ( -z $output ) {
			download_url( $url, $output );
		}
	} else {
		download_url( $url, $output );
	}
	`gzip $output`;
	if ( !( -e $output.".gz" ) ) {
		warn("Couldn't download $output from:\n$url\n");
	}
}

sub get_intergenic_seqs {
	my ( $taxon_oid, $taxon_full_dir ) = @_;
	my $output = $taxon_full_dir.$taxon_oid.".intergenic";

	#return if -e $output."gz";
	my $url = "\"http://img.jgi.doe.gov/cgi-bin/w/main.cgi?section=TaxonDetail&downloadTaxonIntergenicFnaFile=1&_noHeader=1&taxon_oid=".$taxon_oid."\"";
	if ( -e $output.".gz" ) {
		if ( -e $output ) {
			`rm $output`;
		}
		system("gunzip $output.gz");
		if ( -z $output ) {
			download_url( $url, $output );
		}
	} else {
		download_url( $url, $output );
	}
	`gzip $output`;
	if ( !( -e $output.".gz" ) ) {
		warn("Couldn't download $output from:\n$url\n");
	}
}

sub download_url {
	my ( $url, $output ) = @_;
	system("wget -U firefox -q -O $output $url");
	if ( !( -e $output ) ) {
		warn("Couldn't download $output from:\n$url\n");
	} else {
		open( OUT, $output ) || die "Can't open $output for read: $!\n";
		my $head = <OUT>;

		#print $head . "\n";
		if ( $head =~ m/^</ ) {
			warn("Got an html file from:\n$url\n");
		}
		close OUT;
	}
}

1;
