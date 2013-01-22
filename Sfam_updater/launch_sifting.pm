package Sfam_updater::launch_sifting;
use strict;
use Bio::SeqIO;
use Carp;

my $LASTDB_SIZE    = "900M";
my $SLEEP_DURATION = 300;
my $DIV_SIZE       = 1000000;    #how many sequences will be written to 1 file for the sifting;

# Need to launch the Lastal Sifting
# Need to parse lastal Sifting into new family members

sub parse_lastal {
	my %args            = @_;
	my $file_core       = $args{filename_core};
	my $output_dir      = $args{output_dir};
	my @files           = glob( $args{filename_core}."*" );
	my $length_hash_ref = $args{length_hash};
	my %length_hash     = %{$length_hash_ref};
	my %best_hit_hash   = ();
	foreach my $file (@files) {
		print STDERR "parse_lastal Processing $file\n";
		if ( $file =~ m/.gz/ ) {
			## in case of a rerun, need to read from zipped files
			open( IN, "zcat $file |" ) || die "Can't open $file.gz for reading: $!\n";
		} else {
			open( IN, $file ) || die "Can't open $file for reading: $!\n";
		}
		while (<IN>) {
			next if ( $_ =~ /^#/ );    #skipping over comment lines from the lastal output
			next if ( $_ eq "\n" );
			chomp($_);
			my ( $query, $subject, $query_start, $query_end, $bitscore, $frameshift, $subject_start, $subject_end );

			# read table in lastal format
			my @dat = split( /\t/, $_ );
			$bitscore      = $dat[0];
			$subject       = $dat[1];
			$subject_start = $dat[2] + 1;
			$subject_end   = $subject_start + $dat[3] - 1;
			$query         = $dat[6];
			$query_start   = $dat[7] + 1;
			$query_end     = $query_start + $dat[8] - 1;

			#			print "@dat\n";
			next unless defined( $length_hash{$subject} && $length_hash{$query} );
			my $evalue = get_evalue_from_bit_score( length1 => $length_hash{$subject}, length2 => $length_hash{$query}, bitscore => $bitscore );
			my $cov_s  = ( $subject_end - $subject_start ) / $length_hash{$subject};
			my $cov_q  = ( $query_end - $query_start ) / $length_hash{$query};
			if ( $evalue < 1e-10 && $cov_s >= 0.8 && $cov_q >= 0.8 ) {

				#figure out the coverage
				if ( exists $best_hit_hash{$query} && $best_hit_hash{$query}[1] < $bitscore ) {
					@{ $best_hit_hash{$query} } = ();
					push( @{ $best_hit_hash{$query} }, ( $subject, $bitscore ) );
				} else {
					@{ $best_hit_hash{$query} } = ();
					push( @{ $best_hit_hash{$query} }, ( $subject, $bitscore ) );
				}
			}

		}
		close(IN);

		#compress file
		#`gzip $file`;
	}
	open( OUT, ">$output_dir/Lastal_newCDS_to_fam.map" );
	foreach my $new_seq ( keys %best_hit_hash ) {

		#find out what family is the best_hit from
		my $family = Sfam_updater::DB_op::get_family( seq_oid => $best_hit_hash{$new_seq} );

		#insert the new family members
		print OUT "$new_seq\t\$family\n";

		#Sfam_updater::DB_op::insert_new_familymembers(seq_oid => $new_seq, family => $family);
	}
	close(OUT);

	#		print "Added new family members : ".scalar(keys(%best_hit_hash))."\n";

	# Add new familymember

	return "$output_dir/Lastal_newCDS_to_fam.map";
}

sub extract_lastal_leftovers {
	my %args = @_;
	my $family_members_files => $args{family_members_files};
	my $new_cds_core         => $args{new_cds_core};
	my $output_dir           => $args{output_dir};

	#if the output_dir does not exists create it.
	print "Creating $output_dir\n" unless -e $output_dir;
	`mkdir -p $output_dir`         unless -e $output_dir;
	my %new_family_members = ();
	open( FAM_MEM_IN, $family_members_files ) || die "Couldn't open $family_members_files for reading : $!\n";
	while (<FAM_MEM_IN>) {
		next if $_ =~ /#/;    #skip if there is a header;
		$_ =~ m/^(\d+)\s+(\d+)/;
		$new_family_members{$1};
	}
	my $div   = 1;
	my $count = 0;
	open( OUT, ">$output_dir/"."lastal_leftovers_$div.fasta" ) || die "Can't open $output_dir/"."lastal_leftovers_$div.fasta for writing: $!\n";
	my @files = glob( $new_cds_core."*" );
	foreach my $file (@files) {
		print "Processing $file\n";

		#next unless $file =~ m/4/;
		my $inseqs = Bio::SeqIO->new( -file => "$file", -format => 'fasta' );
		while ( my $seq = $inseqs->next_seq() ) {
			my $id       = $seq->display_id();
			my $sequence = $seq->seq();
			if ( !exists $new_family_members{$id} ) {
				$count++;
				if ( $count % $DIV_SIZE == 0 ) {
					close(OUT);
					$div++;
					open( OUT, ">$output_dir/"."lastal_leftovers_$div.fasta" ) || die "Can't open $output_dir/"."lastal_leftovers_$div.fasta for writing: $!\n";
				}
				print OUT ">$id\n$sequence\n";
			}
		}
	}
	return "$output_dir/"."lastal_leftovers";
}

sub get_evalue_from_bit_score {
	my %args     = @_;
	my $length1  = $args{length1};
	my $length2  = $args{length2};
	my $bitscore = $args{bitscore};

	#	print "bitscore : $bitscore\tlength2 : $length2\tlength1 : $length1\t";
	my $evalue = $length1 * $length2 / ( 2**$bitscore );

	#print $evalue."\n";
	my $value = sprintf( "%e", $evalue );

	#print $value."\n";
	return $value;
}

sub compile_sequence_length_hash {
	my %args   = @_;
	my @files1 = glob( $args{file1}."*" );
	my @files2 = glob( $args{file2}."*" );
	print "@files1\n";
	print "@files2\n";
	push( @files1, @files2 );
	my %return_hash = ();
	foreach my $file (@files1) {
		print "Processing $file for length hash\n";

		#next unless $file =~ m/4/;
		my $inseqs = Bio::SeqIO->new( -file => "$file", -format => 'fasta' );
		while ( my $seq = $inseqs->next_seq() ) {
			my $id       = $seq->display_id();
			my $sequence = $seq->seq();
			$return_hash{$id} = length($sequence);
		}
	}
	return \%return_hash;
}

sub index_hmms {
	my %args            = @_;
	my $div_size             = $args{file_size};
	my $username        = $args{username};
	my $password        = $args{password};
	my $db_pointer      = $args{db};
	my $output_dir      = $args{output_dir};
	my $data_repository = $args{repo};
	my $div = 1;
	my $count = 0;
	#if the output_dir does not exists create it.
	print "Creating $output_dir\n" unless -e $output_dir;
	`mkdir -p $output_dir`         unless -e $output_dir;
	my $famids_array_ref = Sfam_updater::DB_op::get_all_famid(
															   db       => $db_pointer,
															   username => $username,
															   password => $password,
	);
	foreach my $fam ( @{$famids_array_ref} ) {
		$count++;
		if ( $count % $div_size == 0 ) {
			close(OUT);
			$div++;
		}
		`zcat $data_repository/$fam.hmm.gz >> $output_dir/HMMs_$div.hmm`;
	}
	return ($count, "$output_dir/HMMs_");
}

sub index_familymembers {
	my %args       = @_;
	my $file       = $args{file_to_index};
	my $output_dir = $args{output_dir};

	#if the output_dir does not exists create it.
	print "Creating $output_dir\n" unless -e $output_dir;
	`mkdir $output_dir`            unless -e $output_dir;

	my $lastdb_command = "lastdb -s $LASTDB_SIZE -p $output_dir/sfams $file";
	my $cmd_output     = `$lastdb_command`;

	if ($cmd_output) {

		#the output of the cmd was not null. Exiting the update.
		croak "Something happened during the indexing. Printing the output\n\n\n $cmd_output\n";
	}
	return "$output_dir/sfams";
}

sub launch_lastal {
	my %args           = @_;
	my $start_index    = 1;
	my $database       = $args{db};
	my $arguments      = $args{arguments};
	my $filename       = $args{filename_core};
	my $output_dir     = $args{output_dir};
	my $error_dir      = $args{error_dir};
	my $lastal_cmd     = "lastal $arguments -o $output_dir/last_sift.ouput.\$SGE_TASK_ID $database $filename"."_\$SGE_TASK_ID.fasta";
 	my @files = glob($filename."_*");
 	my $number_of_jobs = scalar(@files);

	#if the output_dir does not exists create it.
	print "Creating $output_dir\n" unless -e $output_dir;
	`mkdir -p $output_dir`         unless -e $output_dir;

	#if the $error_dir does not exists create it.
	print "Creating $error_dir\n" unless -e $error_dir;
	`mkdir -p $error_dir`         unless -e $error_dir;
	open( OUT, ">$output_dir/lastal_sift.sh" ) || die "Can't open $output_dir/lastal_sift.sh for writing: $!\n";
	print OUT "
#!/bin/sh
#\$ -cwd
#\$ -V
#\$ -S /bin/bash

#\$ -t $start_index-$number_of_jobs
#\$ -e $error_dir/lastal_sifting_jobs.err
#\$ -o $error_dir/lastal_sifting_jobs.out

$lastal_cmd 

";
	close(OUT);
	`chmod 755 $output_dir/lastal_sift.sh`;
	my $qsub_command = "qsub -q eisen.q $output_dir/lastal_sift.sh";
	my $qsub         = `$qsub_command`;
	$qsub =~ /Your job (\d+) /;
	my $job_id = $1;
	my $flag   = 1;

	while ($flag) {
		my $output = `qstat -j $job_id 2>&1`;
		if ( !defined($output) || $output =~ /Following jobs do not exist/ ) {
			$flag = 0;
		}
		sleep($SLEEP_DURATION);
	}
	return "$output_dir/last_sift.ouput";
}
