package Sfam_updater::launch_sifting;
use strict;

#use Sfam_updater::DB_op;
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
		`gzip $file`;
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
	my $output_file          => $args{output_file};

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
	my $count = 0;
	open( OUT, ">$output_dir/"."lastal_leftovers.fasta" ) || die "Can't open $output_dir/"."lastal_leftovers.fasta for writing: $!\n";
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
				print OUT ">$id\n$sequence\n";
			}
		}
	}
	return "$output_dir/"."lastal_leftovers";
}

sub extract_hmmsearch_leftovers {
	my %args                = @_;
	my $family_members_file = $args{family_members_files};
	my $new_cds_file        = $args{new_cds_file};
	my $output_dir          = $args{output_dir};

	#if the output_dir does not exists create it.
	print "Creating $output_dir\n" unless -e $output_dir;
	`mkdir -p $output_dir`         unless -e $output_dir;
	my %family_members = ();
	open( IN_FAM_MAP, $family_members_file ) || die "Can't open $family_members_file for reading: $!\n";
	while (<IN_FAM_MAP>) {
		next if $_ =~ m/^#/;
		$_ =~ m/^(\S+)\s+(\S+)/;
		$family_members{$1} = 1;
	}
	close(IN_FAM_MAP);
	my $count = 0;
	open( OUT, ">$output_dir/hmmsearch_leftovers.fasta" ) || die "Can't open $output_dir/hmmsearch_leftovers.fasta for reading: $!\n";
	my $inseqs = Bio::SeqIO->new( -file => "$new_cds_file", -format => 'fasta' );
	while ( my $seq = $inseqs->next_seq() ) {
		my $id       = $seq->display_id();
		my $sequence = $seq->seq();
		if ( !exists $family_members{$id} ) {
			$count++;
			print OUT ">$id\n$sequence\n";
		}
	}
	close(OUT);
	print STDERR "found $count left over sequences after hmmsearch\n";

	#split the hmmsearch leftovers into smaller chunks
	my $div_count = int( $count / 14 );
	$count = 0;
	my $div = 1;
	open( OUT, ">$output_dir/hmmsearch_leftovers_$div.fasta" ) || die "Can't open $output_dir/hmmsearch_leftovers_$div.fasta for writing: $!\n";
	my $inseqs = Bio::SeqIO->new( -file => "$output_dir/hmmsearch_leftovers.fasta", -format => 'fasta' );
	while ( my $seq = $inseqs->next_seq() ) {
		$count++;
		my $id       = $seq->display_id();
		my $sequence = $seq->seq();
		if ( $count % 14 ) {
			close(OUT);
			$div++;
			open(">$output_dir/hmmsearch_leftovers_$div.fasta") || die "Can't open $output_dir/hmmsearch_leftovers_$div.fasta for writing: $!\n";
		}
		print OUT ">$id\n$sequence\n";
	}
	return $output_dir."/hmmsearch_leftovers";
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

sub parse_hmmsearch {
	my %args       = @_;
	my $file_core  = $args{filename_core};
	my $output_dir = $args{output_dir};
	my @files      = glob( $args{filename_core}.".*" );
	my %top_both   = ();
	print "$output_dir\n$args{filename_core}\n\n@files\n";
	foreach my $file (@files) {
		print STDERR "parse_hmmsearch Processing $file\n";
		if ( $file =~ m/.gz/ ) {
			## in case of a rerun, need to read from zipped files
			open( IN, "zcat $file |" ) || die "Can't open $file.gz for reading: $!\n";
		} else {
			open( IN, $file ) || die "Can't open $file for reading: $!\n";
		}
		while (<IN>) {
			next if ( $_ =~ /^#/ );    #skipping over comment lines from the lastal output
			chomp($_);
			my ( $query, $subject, $query_start, $query_end, $bitscore, $frameshift, $subject_start, $subject_end );
			next unless $_ =~ m/(\S+)\s+-\s+(\d+)\s+(\d+)\s+-\s+(\d+)\s+(\S+)\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)/;
			my $seq      = $1;
			my $seqLen   = $2;
			my $hmm      = $3;
			my $hmmLen   = $4;
			my $evalue   = $5;
			my $hmmStart = $6;
			my $hmmEnd   = $7;
			my $seqStart = $8;
			my $seqEnd   = $9;
			my $covSeq   = ( $seqEnd - $seqStart ) / $seqLen;
			my $covHmm   = ( $hmmEnd - $hmmStart ) / $hmmLen;

			#           print $hmm."\n";
			#$hmm =~ s/[^\d]//g;
			#           print $hmm."\n";
			#they are both satisfying the coverage threshold
			if ( $covSeq >= 0.8 && $covHmm >= 0.8 ) {
				if ( exists $top_both{$seq} ) {
					if ( $evalue < $top_both{$seq}[1] ) {
						delete( $top_both{seq} );
						push( @{ $top_both{$seq} }, ( $hmm, $evalue ) );
					} else {

						#do nothing
					}
				} else {
					push( @{ $top_both{$seq} }, ( $hmm, $evalue ) );
				}
			}
		}
		close(IN);

		#compress file
		`gzip $file` if $file !~ m/.gz/;
	}
	open( OUT, ">$output_dir/hmmsearch_newCDS_to_fam.map" );
	foreach my $new_seq ( keys %top_both ) {

		#insert the new family members
		print OUT "$new_seq\t\$top_both{$new_seq}[0]\n";

		#Sfam_updater::DB_op::insert_new_familymembers(seq_oid => $new_seq, family => $family);
	}
	close(OUT);

	return "$output_dir/hmmsearch_newCDS_to_fam.map";
}

sub compile_sequence_length_hash {
	my %args   = @_;
	my @files1 = glob( $args{file1}."*" ) if defined $args{file1};
	my @files2 = glob( $args{file2}."*" ) if defined $args{file2};
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
	my $div_size        = $args{file_size};
	my $number_of_files = $args{number_files};
	my $username        = $args{username};
	my $password        = $args{password};
	my $db_pointer      = $args{db};
	my $output_dir      = $args{output_dir};
	my $data_repository = $args{repo};
	my $div             = 1;
	my $count           = 0;

	#if the output_dir does not exists create it.
	print "Creating $output_dir\n" unless -e $output_dir;
	`mkdir -p $output_dir`         unless -e $output_dir;
	my $famids_array_ref = Sfam_updater::DB_op::get_all_famid(
															   db       => $db_pointer,
															   username => $username,
															   password => $password,
	);
	$div_size = int( scalar @{$famids_array_ref} / $number_of_files ) if defined $number_of_files;
	foreach my $fam ( @{$famids_array_ref} ) {
		my @array_fam = @{$fam};
		$count++;
		if ( $count % $div_size == 0 ) {
			close(OUT);
			$div++;
		}
		print STDERR "zcat $data_repository/$array_fam[0].hmm.gz >> $output_dir/HMMs_$div.hmm\n";
		`zcat $data_repository/$array_fam[0].hmm.gz >> $output_dir/HMMs_$div.hmm`;
	}
	return ( $count, $output_dir."/HMMs" );
}

#need to Get new famID and use that when writing out the seqID to famID mapping
sub parse_mcl {
	my %args            = @_;
	my $mcl_output_file = $args{mcl_output_file};
	my $output_dir      = $args{output_dir};
	my $db              = $args{db};
	my $username        = $args{username};
	my $password        = $args{password};
	my $output_dir      = $args{output_dir};
	print "Creating $output_dir\n" unless -e $output_dir;
	`mkdir -p $output_dir`         unless -e $output_dir;
	my $analysis = MRC->new();
	$analysis->set_dbi_connection($db);
	$analysis->set_username($username);
	$analysis->set_password($password);
	$analysis->build_schema();
	open( MCL_IN, $mcl_output_file ) || die "Can't open $mcl_output_file for reading: $!\n";
	open( FAM_MAPPING, ">$output_dir/mcl_newCDS_to_fam.map" ) || "Can't open $output_dir/mcl_newCDS_to_fam.map for writing: $!\n";
	my $familyconstruction_id = $analysis->MRC::DB::get_max_familyconstruction_id();
	print STDERR "Family construction ID used : $familyconstruction_id\n";

	while (<MCL_IN>) {
		chomp($_);
		my @gene_oids = split( /\t/, $_ );      #1 family per line, gene_oids separated by tabs
		my $family_size = scalar(@gene_oids);
		next if $family_size <= 1;              #ignore singletons
		my (
			 $fam_alt_id, $name,          $description,      $alnpath,  $seed_alnpath, $hmmpath,
			 $reftree,    $alltree,       $universality,     $evenness, $arch_univ,    $bact_univ,
			 $euk_univ,   $unknown_genes, $pathogen_percent, $aquatic_percent
		  )
		  = undef;
		$analysis->MRC::DB::insert_family(
										   $familyconstruction_id, $fam_alt_id,       $name,      $description, $alnpath,
										   $seed_alnpath,          $hmmpath,          $reftree,   $alltree,     $family_size,
										   $universality,          $evenness,         $arch_univ, $bact_univ,   $euk_univ,
										   $unknown_genes,         $pathogen_percent, $aquatic_percent
		);
		my $family_ID = $analysis->MRC::DB::get_max_famid();
		foreach my $seqID (@gene_oids) {
			print FAM_MAPPING $seqID."\t".$family_ID."\n";
		}
	}
	close(MCL_IN);
	close(FAM_MAPPING);

	#insert into family
	#get max famID
	#use famiID to write out familymembers to file
	# return mapping filename
	return "$output_dir/mcl_newCDS_to_fam.map";
}

sub fine_tune_representative_picking {
	my %args      = @_;
	my $blast_dir = $args{blast_dir};
	my @files     = <$blast_dir/*.abc>;
	foreach my $file (@files) {
		next if $file !~ m/\/(\d+).abc/;
		my $core = $1;
		if ( !-z $file ) {
			if ( !-e "$blast_dir/$core.mcl.log" ) {
				print "Initializing reps for $core\n";
				`perl Sfam_updater/pick_rep_by_mcl.pl -i $blast_dir/$core.abc -o $blast_dir/$core.mcl -c 99`;
			}
			print "Done initializing reps $core\n";
			my $rep_count = 10000;
			while ( $rep_count > 250 ) {
				open( inFILE, "$blast_dir/$core.mcl.log" );
				while (<inFILE>) {
					if ( $_ =~ m/number_of_rep/ ) {

						#do nothing
					} elsif ( $_ =~ m/Final: (\d+)\s+(\d+)/ ) {
						$rep_count = $1;
						my $cutoff = $2;
						if ( $rep_count > 1500 ) {
							$cutoff = $cutoff - 5;
						} elsif ( $rep_count < 250 ) {
							print "Reps # under the cutoff, skipping\n";
							next;
						} else {
							$cutoff = $cutoff - 1;
						}
						`perl Sfam_updater/pick_rep_by_mcl.pl -i $blast_dir/$core.parsed -o $blast_dir/$core.mcl -c 99`;
					}
				}
			}
		}
	}
}

sub launch_mcl {
	my %args         = @_;
	my $queue_params = $args{queue};
	my $output_dir   = $args{output_dir};
	my $input_file   = $args{input};
	my $mcl_params   = $args{mcl_params};
	my $error_dir    = $args{error_dir};
	my $threads      = $args{threads};

	#if the output_dir does not exists create it.
	print "Creating $output_dir\n" unless -e $output_dir;
	`mkdir -p $output_dir`         unless -e $output_dir;
	my $mcl_cmd = "mcl $input_file $mcl_params -te $threads --abc -o $output_dir/mcl_ouput.mcl";
	open( OUT, ">$output_dir/mcl_job.sh" ) || die "Can't open $output_dir/mcl_job.sh for writing: $!\n";
	print OUT "
#!/bin/sh
#\$ -cwd
#\$ -V
#\$ -S /bin/bash

#\$ -e $error_dir/mcl_job.err
#\$ -o $error_dir/mcl_job.out

$mcl_cmd

";
	close(OUT);
	`chmod 755 $output_dir/mcl_job.sh`;
	my $qsub_command = "qsub $queue_params -pe threaded $threads $output_dir/mcl_job.sh";
	my $qsub         = `$qsub_command`;
	$qsub =~ /Your job (\d+) /;
	my $job_id = $1;
	my $flag   = 1;

	while ($flag) {
		sleep($SLEEP_DURATION);
		my $output = `qstat -j $job_id 2>&1`;
		if ( !defined($output) || $output =~ /Following jobs do not exist/ ) {
			$flag = 0;
		}
	}
	return "$output_dir/mcl_ouput.mcl";
}

sub parse_blast {
	my %args               = @_;
	my $files_to_parse     = $args{blast_results_core};
	my $coverage           = $args{coverage};
	my $evalue             = $args{evalue};
	my $output_dir         = $args{output_dir};
	my $seq_length_hashref = $args{seq_lengths};
	my %seq_lengths        = %{$seq_length_hashref};

	#if the output_dir does not exists create it.
	print "Creating $output_dir\n" unless -e $output_dir;
	`mkdir -p $output_dir`         unless -e $output_dir;
	my @files = glob( $files_to_parse."*" );
	open( OUT, ">$output_dir/mcl_input.abc" ) || die "Can't open $output_dir/mcl_input.abc for writing: $!\n";
	foreach my $file (@files) {
		open( IN, $file ) || die "Can't open $file for reading: $!\n";
		while (<IN>) {
			next if $_ =~ m/^#/;
			my @line = split( /\t/, $_ );
			next if $line[0] eq $line[1];
			my $queryID   = $line[0];
			my $hitID     = $line[1];
			my $hit_len   = $seq_lengths{$hitID};
			my $query_len = $seq_lengths{$queryID};
			next unless defined($hit_len);
			next unless defined($query_len);
			my $eValue     = $line[10];
			my $hitStart   = $line[8];
			my $hitEnd     = $line[9];
			my $queryStart = $line[6];
			my $queryEnd   = $line[7];
			my $identity   = $line[2] / 100;
			$identity = sprintf( "%.3f", $identity );
			my $query_coverage = ( $queryEnd - $queryStart + 1 ) / $query_len;
			my $hit_coverage   = ( $hitEnd - $hitStart + 1 ) / $hit_len;

			#	print "$query_coverage\t$hit_coverage\t$eValue\n";

			if (    ( $query_coverage > $coverage )
				 && ( $hit_coverage > $coverage )
				 && ( $eValue < $evalue ) )
			{
				print OUT "$queryID\t$hitID\t$identity\n";
			}
		}
		close(IN);
	}
	close(OUT);
	return "$output_dir/mcl_input.abc";
}

sub launch_blast {
	my %args                   = @_;
	my $blast_input_files_core = $args{blast_input_files_core};
	my $arguments              = $args{arguments};
	my $output_dir             = $args{output_dir};
	my $error_dir              = $args{error_dir};
	my $threads                = $args{threads};
	my %job_ids                = ();

	#if the output_dir does not exists create it.
	print "Creating $output_dir\n" unless -e $output_dir;
	`mkdir -p $output_dir`         unless -e $output_dir;

	#if the $error_dir does not exists create it.
	print "Creating $error_dir\n" unless -e $error_dir;
	`mkdir -p $error_dir`         unless -e $error_dir;
	my @files       = glob( $blast_input_files_core."*" );
	my $file_number = scalar(@files);
	print STDERR "$blast_input_files_core\n";
	print STDERR "There are $file_number Blast files\n";
	print "ARGUMENTS $arguments\n";
	my $count = 0;

	for ( my $i = 1; $i <= $file_number; $i++ ) {
		$count++;
		my $blast_cmd = "blastp -num_threads $threads $arguments ";
		$blast_cmd .= "-subject $blast_input_files_core"."_$i.fasta ";
		$blast_cmd .= "-query $blast_input_files_core"."_\$SGE_TASK_ID.fasta ";
		$blast_cmd .= "-out $output_dir/blast_output.$i.\$SGE_TASK_ID.tblout";

		open( OUT, ">$output_dir/all_v_all_blast_$count.sh" ) || die "Can't open $output_dir/all_v_all_blast_$count.sh: $!\n";
		print OUT "
#!/bin/sh
#\$ -cwd
#\$ -V
#\$ -S /bin/bash

#\$ -t 1-$file_number
#\$ -e $error_dir/all_v_all_jobs.err
#\$ -o $error_dir/all_v_all_jobs.out

$blast_cmd

";
		close(OUT);
		`chmod 755 $output_dir/all_v_all_blast_$count.sh`;
		my $qsub_command = "qsub -q eisen.q -pe threaded $threads $output_dir/all_v_all_blast_$count.sh";
		my $qsub         = `$qsub_command`;
		$qsub =~ /Your job (\d+) /;
		$job_ids{$1} = 1;

	}
	my $flag = 1;

	while ( scalar( keys(%job_ids) ) > 1 ) {
		foreach my $jobid ( keys(%job_ids) ) {
			my $output = `qstat -j $jobid 2>&1`;
			if ( !defined($output) || $output =~ /Following jobs do not exist/ ) {
				delete $job_ids{$jobid};
			}
		}
		sleep($SLEEP_DURATION);
	}
	return "$output_dir/hmmsearch_sift.ouput";
}

sub launch_hmmsearch {
	my %args              = @_;
	my $hmm_files_core    = $args{hmmfiles};
	my $seqs              = $args{seq_files};
	my $arguments         = $args{arguments};
	my $output_dir        = $args{output_dir};
	my $error_dir         = $args{error_dir};
	my $hmmsearch_threads = $args{threads};

	#if the output_dir does not exists create it.
	print "Creating $output_dir\n" unless -e $output_dir;
	`mkdir -p $output_dir`         unless -e $output_dir;

	#if the $error_dir does not exists create it.
	print "Creating $error_dir\n" unless -e $error_dir;
	`mkdir -p $error_dir`         unless -e $error_dir;
	my @files       = glob( $hmm_files_core."*" );
	my $file_number = scalar(@files);
	print STDERR "There are $file_number HMMs\n";
	$file_number = 3;
	print "ARGUMENTS $arguments\n";

	my $hmmsearch_cmd =
	  "hmmsearch --cpu $hmmsearch_threads $arguments --domtblout $output_dir/hmmsearch_sift.ouput.\$SGE_TASK_ID $hmm_files_core"."_\$SGE_TASK_ID.hmm $seqs";
	open( OUT, ">$output_dir/hmmsearch_sift.sh" ) || die "Can't open $output_dir/lastal_sift.sh for writing: $!\n";
	print OUT "
#!/bin/sh
#\$ -cwd
#\$ -V
#\$ -S /bin/bash

#\$ -t 1-$file_number
#\$ -e $error_dir/hmmsearch_sifting_jobs.err
#\$ -o $error_dir/hmmsearch_sifting_jobs.out

$hmmsearch_cmd

";
	close(OUT);
	`chmod 755 $output_dir/hmmsearch_sift.sh`;
	my $qsub_command = "qsub -q eisen.q -pe threaded $hmmsearch_threads $output_dir/hmmsearch_sift.sh";
	my $qsub         = `$qsub_command`;
	$qsub =~ /Your job (\d+) /;
	my $job_id = $1;
	my $flag   = 1;

	while ($flag) {
		sleep($SLEEP_DURATION);
		my $output = `qstat -j $job_id 2>&1`;
		if ( !defined($output) || $output =~ /Following jobs do not exist/ ) {
			$flag = 0;
		}
	}
	return "$output_dir/hmmsearch_sift.ouput";
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
	my @files          = glob( $filename."_*" );
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
		sleep($SLEEP_DURATION);
		my $output = `qstat -j $job_id 2>&1`;
		if ( !defined($output) || $output =~ /Following jobs do not exist/ ) {
			$flag = 0;
		}
	}
	return "$output_dir/last_sift.ouput";
}

sub build_aln_hmm_trees {
	my %args            = @_;
	my $type            = $args{type};
	my $directory       = $args{directory};
	my $output_dir      = $args{output};
	my $error_dir       = $args{error};
	my $total_jobs      = $args{total_jobs};
	my $data_repository = $args{repo};
	##if the output_dir does not exists create it.
	print "Creating $output_dir\n" unless -e $output_dir;
	`mkdir -p $output_dir`         unless -e $output_dir;
	print $directory."/*_$type"."CDS.fasta\n";
	my @files = glob( $directory."/*_$type"."CDS.fasta" );
	##if the error_dir does not exists create it.
	print "Creating $error_dir\n" unless -e $error_dir;
	`mkdir -p $error_dir`         unless -e $error_dir;
	my $low_file  = undef;
	my $high_file = 0;
	foreach my $file (@files) {
		$file =~ m/\/(\d+)_.*CDS.fasta/;
		$low_file  = $1 if ( !defined $low_file );
		$low_file  = $1 if $low_file > $1;
		$high_file = $1 if $high_file < $1;
	}
	my $end_array_job = $low_file + $total_jobs - 1;
	my $cmd = "perl external_software_launcher.pl $type \$SGE_TASK_ID $high_file $total_jobs $directory _$type"."CDS.fasta $output_dir $data_repository";
	print "$cmd\n";
	print "Lowfile :$low_file\tHigh file : $high_file\n";
	open( OUT, ">$output_dir/build_aln_hmm_trees.sh" ) || die "Can't open $output_dir/build_aln_hmm_trees.sh for writing: $!\n";
	print OUT "
#!/bin/sh
#\$ -cwd
#\$ -V
#\$ -S /bin/bash

#\$ -t $low_file-$end_array_job
#\$ -e $error_dir/lastal_sifting_jobs.err
#\$ -o $error_dir/lastal_sifting_jobs.out

$cmd 

";
	close(OUT);
	print "RUNNING $cmd\n";
	`chmod 755 $output_dir/build_aln_hmm_trees.sh`;
	my $qsub_command = "qsub -q eisen.q $output_dir/build_aln_hmm_trees.sh";
	my $qsub         = `$qsub_command`;
	$qsub =~ /Your job (\d+) /;
	my $job_id = $1;
	my $flag   = 1;

	while ($flag) {
		sleep($SLEEP_DURATION);
		my $output = `qstat -j $job_id 2>&1`;
		if ( !defined($output) || $output =~ /Following jobs do not exist/ ) {
			$flag = 0;
		}
	}
	return 1;

	#for(my $i = 1; $i <= $step; $i++){
	#}
}

