#! /usr/bin/perl -w
use strict;
use File::Basename;
use Getopt::Long;

use GTF;
use EGlite;

# hashes I will use later to ambiguate and disambiguate pairs of base calls
my %ambig_hash;
populate_ambig_hash(\%ambig_hash);
my %disambig_hash;
populate_disambig_hash(\%disambig_hash);
my %type_hash;
populate_type_hash(\%type_hash);

my $mode;
GetOptions( "mode=s" => \$mode);

# input files
my ($trace_file, $top_sequence_file) = @ARGV;

my $num_of_arguments = @ARGV;

my $usage = "trace_recalling.pl [--mode=single,<numeric threshold>] <trace file> <genomic sequence file>";

if($num_of_arguments != 2) {
    die($usage);
}

my ($trace_base, $trace_dir, $trace_ext) = fileparse($trace_file,'\..*');

my $trace_base_and_ext = $trace_base.$trace_ext;

my $expt_dir = $trace_dir.$trace_base_and_ext."_dir";

if(!(-e $expt_dir)) {
    system("mkdir $expt_dir");
}

my $basename = $expt_dir."/".$trace_base.$trace_ext;

if(defined($mode)) {
    
    if($mode =~ /single,([\d\.]+)/) {
	my $threshold = $1;

	my $sequence_file = $top_sequence_file;


        # build neccessary filenames for later
	my $poly                = $expt_dir."/".$trace_base_and_ext.".poly";
	my $ambig_file          = $basename.".ambig";
	my $recall_file         = $basename.".recall";
	my $first_align         = $basename.".first_align";
	my $second_align        = $basename.".second_align";	


	my $reverse_complement_flag;
	
# get the original poly file and make an internal representation of it
	system( "phred -dd $expt_dir $trace_file");
	
# create poly hash and set initial threshold for calling ambiguity codes
	my %poly_hash;
	create_poly_hash($poly, \%poly_hash, 1/$threshold);
	
# generate ambiguity sequence the way original trace recalling does it (no thresholding at all so far)
	my $ambig_seq = blindly_generate_ambig_sequence(\%poly_hash, \%ambig_hash);
	
# dump ambiguity sequence to a file
	open(AMBIG, ">$ambig_file");
	print AMBIG ">${ambig_file}\n";
	print AMBIG $ambig_seq, "\n";   
	close(AMBIG);
	
# run est_genome_ambig to generate the first align file
	system( "est_genome_ambig -genome $sequence_file -est $ambig_file -mismatch 6 -gap_panalty 6 -match 2 -splice_penalty 20 -intron_penalty 40 -align -minscore 0 > $first_align");
	
	my $first_align_report = get_estgenome_report($first_align);
	
	if($first_align_report->genomeDirection eq "forward" && $first_align_report->estDirection eq "forward") {
	    $reverse_complement_flag = 0;
	}
	elsif($first_align_report->genomeDirection eq "forward" && $first_align_report->estDirection eq "reversed") {
	    
	    $reverse_complement_flag = 1;
	    
	    # in this case need to reverse complement the sequence to make alignment foward to foward
	    my $new_sequence_file = $sequence_file."_rc";
	    system("reverse_complement.pl $sequence_file > $new_sequence_file");
	    $sequence_file = $new_sequence_file;
	    
	    system( "est_genome_ambig -mode forward -genome $sequence_file -est $ambig_file -mismatch 6 -gap_panalty 6 -match 2 -splice_penalty 20 -intron_penalty 40 -align -minscore 0 > $first_align");
	    
	}
	else {
	    die("with first alignment... doesnt appear to be fwd/fwd or rev/fwd and should be one of those 2!");
	}
	
	# parse first alignment into a hash
	my %first_align_contents;
	create_alignment_hash($first_align, \%first_align_contents);
	
	#first_align_to_trace takes a position in the first alignment and returns the corresponding position in the abmiguity sequence    
	my %first_align_to_trace;
	
	# trace_to_first_align takes a postion in the ambiguity sequence and returns the corresponding position in the first alignment
	my %trace_to_first_align;
	
	# create mapping between the trace (ambiguity sequence) and first alignment... 
	create_sequence_alignment_map(\%trace_to_first_align, \%first_align_to_trace, \%first_align_contents, $poly_hash{"length"});
	
	# generate the recalled sequence from the first align file
	my $recall_seq = blindly_generate_recall_sequence(\%poly_hash, \%first_align_contents, \%trace_to_first_align, \%disambig_hash);
	
	# dump recall sequence out to a file
	open(RECALL, ">$recall_file");
	print RECALL ">${recall_file}\n";
	print RECALL $recall_seq, "\n";   
	close(RECALL);
	
	system( "est_genome_ambig -genome $sequence_file -est $recall_file -match 1 -mismatch 1 -gap_panalty 2 -splice_penalty 20 -intron_penalty 40 -align -minscore 0 > $second_align");	
		
	exit(0);
    }
    else {
	die("Inivalid mode");
    }
}

# if you get to this point it means go ahead and run the analysis as in the paper

my @logfile_list;
my $movie_log = $basename.".result_log";

# default here should be up to 20
# the whole program basically takes place in this loop
for(my $thresh = 1; $thresh <= 20; $thresh++) {

    print STDERR "working on threshold $thresh\n";

    my $sequence_file = $top_sequence_file;

# build neccessary filenames for later

    my $poly = $expt_dir."/".$trace_base_and_ext.".poly";

    my $ambig_file          = $basename.".ambig".".${thresh}";
    my $recall_file         = $basename.".recall".".${thresh}";
    my $first_align         = $basename.".first_align".".${thresh}";
    my $first_gtf           = $basename.".first_gtf".".${thresh}";
    my $second_align        = $basename.".second_align".".${thresh}";
    my $second_gtf          = $basename.".second_gtf".".${thresh}";

    my $ambig_file2         = $basename.".ambig2".".${thresh}";
    my $recall_file2        = $basename.".recall2".".${thresh}";
    my $first_align2        = $basename.".first_align2".".${thresh}";
    my $first_gtf2          = $basename.".first_gtf2".".${thresh}";
    my $second_align2       = $basename.".second_align2".".${thresh}";
    my $second_gtf2         = $basename.".second_gtf2".".${thresh}";
    my $first_pass_log_file = $basename.".first_log".".${thresh}";
    my $log_file            = $basename.".iter_log".".${thresh}";
    
    my $acedb_command_line  = $basename.".acedb_command".".${thresh}";

    my $reverse_complement_flag;
    
    push(@logfile_list, $log_file);
    
# get the original poly file and make an internal representation of it
    system( "phred -dd $expt_dir $trace_file");
    
# create poly hash and set initial threshold for calling ambiguity codes
    my %poly_hash;
    create_poly_hash($poly, \%poly_hash, $thresh);
    
# generate ambiguity sequence the way original trace recalling does it (no thresholding at all so far)
    my $ambig_seq = blindly_generate_ambig_sequence(\%poly_hash, \%ambig_hash);
    
# dump ambiguity sequence to a file
    open(AMBIG, ">$ambig_file");
    print AMBIG ">${ambig_file}\n";
    print AMBIG $ambig_seq, "\n";   
    close(AMBIG);
    
# run est_genome_ambig to generate the first align file
    system( "est_genome_ambig -genome $sequence_file -est $ambig_file -mismatch 6 -gap_panalty 6 -match 2 -splice_penalty 20 -intron_penalty 40 -align -minscore 0 > $first_align");
    
    my $first_align_report = get_estgenome_report($first_align);
    
    if($first_align_report->genomeDirection eq "forward" && $first_align_report->estDirection eq "forward") {
	$reverse_complement_flag = 0;
    }
    elsif($first_align_report->genomeDirection eq "forward" && $first_align_report->estDirection eq "reversed") {

	$reverse_complement_flag = 1;

        # in this case need to reverse complement the sequence to make alignment foward to foward
	my $new_sequence_file = $sequence_file."_rc";
	system("reverse_complement.pl $sequence_file > $new_sequence_file");
	$sequence_file = $new_sequence_file;

	system( "est_genome_ambig -mode forward -genome $sequence_file -est $ambig_file -mismatch 6 -gap_panalty 6 -match 2 -splice_penalty 20 -intron_penalty 40 -align -minscore 0 > $first_align");
	
    }
    else {
	die("with first alignment... doesnt appear to be fwd/fwd or rev/fwd and should be one of those 2!");
    }

    # convert first align file to gtf    
    system("est2geno_to_gtf.pl $first_align > $first_gtf");
    
    # parse first alignment into a hash
    my %first_align_contents;
    create_alignment_hash($first_align, \%first_align_contents);

    #first_align_to_trace takes a position in the first alignment and returns the corresponding position in the abmiguity sequence    
    my %first_align_to_trace;

    # trace_to_first_align takes a postion in the ambiguity sequence and returns the corresponding position in the first alignment
    my %trace_to_first_align;
    
    # create mapping between the trace (ambiguity sequence) and first alignment... 
    create_sequence_alignment_map(\%trace_to_first_align, \%first_align_to_trace, \%first_align_contents, $poly_hash{"length"});
    
    # generate the recalled sequence from the first align file
    my $recall_seq = blindly_generate_recall_sequence(\%poly_hash, \%first_align_contents, \%trace_to_first_align, \%disambig_hash);
    
    # dump recall sequence out to a file
    open(RECALL, ">$recall_file");
    print RECALL ">${recall_file}\n";
    print RECALL $recall_seq, "\n";   
    close(RECALL);
    
    # run est_genome on the recalled sequence and the genomic sequence
    system( "est_genome_ambig -genome $sequence_file -est $recall_file -match 1 -mismatch 1 -gap_panalty 2 -splice_penalty 20 -intron_penalty 40 -align -minscore 0 > $second_align");	
 
    my $second_align_report = get_estgenome_report($second_align);
    
    if($second_align_report->genomeDirection ne "forward" || $second_align_report->estDirection ne "forward") {
	open(LOG, ">${log_file}");
	close(LOG);
	next;
    }
    

    system("est2geno_to_gtf.pl $second_align > $second_gtf");
    
    open(LOG, ">$first_pass_log_file");
    
    my $first_compare_string = `compare_altsplices_single.pl ${first_gtf} ${second_gtf} --no_both_altsplice`;
    
    if($first_compare_string ne "") {
	print LOG "altsplice found --> $first_compare_string\n";
    }
    else {
	print LOG "no altsplice found\n";
    }
    
    close(LOG);
    
    # parse second alignment into a has in the same way as the first alignment
    my %second_align_contents;
    create_alignment_hash($second_align, \%second_align_contents);
    
    my %second_align_to_recall;
    my %recall_to_second_align;

    create_sequence_alignment_map(\%recall_to_second_align, \%second_align_to_recall, \%second_align_contents, $poly_hash{"length"});

    # walk through the poly file and take a look at what happened at each position..
    tester(\%poly_hash, \%trace_to_first_align, \%first_align_contents, \%recall_to_second_align, \%second_align_contents, \%type_hash, $thresh);
    
    my($new_ambig, $new_recall) = guided_generate_sequences(\%poly_hash);
    
    
    # dump new ambiguity sequence to a file
    open(AMBIG, ">$ambig_file2");
    print AMBIG ">${ambig_file2}\n";
    print AMBIG $new_ambig, "\n";   
    close(AMBIG);
    
    # dump new recall sequence out to a file
    open(RECALL, ">$recall_file2");
    print RECALL ">${recall_file2}\n";
    print RECALL $new_recall, "\n";   
    close(RECALL);
    
# run est_genome_ambig to generate the first align file
# using this one

    system( "est_genome_ambig -genome $sequence_file -est $ambig_file2 -mismatch 6 -gap_panalty 6 -match 2 -splice_penalty 20 -intron_penalty 40 -align -minscore 0 > $first_align2");

    system("est2geno_to_gtf.pl $first_align2 > $first_gtf2");
    
    # run est_genome on the recalled sequence and the genomic sequence... this is the old secondary alignment 

    system( "est_genome_ambig -genome $sequence_file -est $recall_file2 -match 1 -mismatch 1 -gap_panalty 4 -splice_penalty 20 -intron_penalty 40 -align -minscore 0 > $second_align2");

    system("est2geno_to_gtf.pl $second_align2 > $second_gtf2");
    
    # once you have the first and second second pass gtf files check them
    # for alternate splices and dump the results into a log file
    
    open(LOG, ">$log_file");
    
    my $compare_string = `compare_altsplices_single.pl ${first_gtf2} ${second_gtf2} --no_both_altsplice`;
    
    if($compare_string ne "") {
	print LOG "altsplice found --> $compare_string\n";
    }
    else {
	print LOG "no altsplice found\n";
    }
    
    close(LOG);
}


# look at log files and alignments of different thresholds and try and figure out if you really have an alternate splice
analyze_movie(\@logfile_list, $movie_log);

exit(0);

################## SUBROUTINES ##############################

sub make_add_string {

    my($gtf_file1, $gtf_file2) = @_;

    my $gtf1 = GTF::new({gtf_filename => $gtf_file1});
    my $genes1 = $gtf1->transcripts;
    
    my $gtf2 = GTF::new({gtf_filename => $gtf_file2});
    my $genes2 = $gtf2->transcripts;
    
    my $max = 0;
    my $min = 1000000000000000;
    
    my $length1 = 0;
    my $length2 = 0;
    
    foreach my $gene (@$genes1) {
	
	my $cds = $gene->cds;
	foreach my $c (@$cds){
	    
	    $length1 += $c->stop - $c->start;
	    
	    if($c->stop > $max) {
		$max = $c->stop;
	    }
	    if($c->start < $min) {
		$min = $c->start;
	    }
	    
	}
	
    }
    
    foreach my $gene (@$genes2){
	
	my $cds = $gene->cds;
	foreach my $c (@$cds){
	    $length2 += $c->stop - $c->start;
	    
	    
	    if($c->stop > $max) {
		$max = $c->stop;
	    }
	    if($c->start < $min) {
		$min = $c->start;
	    }
	}
    }
    
    
    my $compare_string1 = "";
    my $compare_string2 = "";
    my $add_string = "";
    
    for(my $i=0;$i<$max - $min;$i++) {
	
	$compare_string1 .= "0";
	$compare_string2 .= "0";
    }
    
    
    foreach my $gene (@$genes1){
	
	my $cds = $gene->cds;
	foreach my $c (@$cds){
	    
	    my $temp = "";
	    for(my $i=0;$i<$c->stop - $c->start;$i++) {
		$temp .= "1";
	    }
	    substr($compare_string1, $c->start - $min, $c->stop - $c->start, $temp);
	}
	
    }
    
    foreach my $gene (@$genes2){
	
	my $cds = $gene->cds;
	foreach my $c (@$cds){
	    
	    my $temp = "";
	    for(my $i=0;$i<$c->stop - $c->start;$i++) {
		$temp .= "1";
	    }
	    substr($compare_string2, $c->start - $min, $c->stop - $c->start, $temp);
	    
	}
    }
    
    for(my $i=0;$i<length($compare_string1);$i++) {
	
	my $c1 = substr($compare_string1, $i, 1);
	my $c2 = substr($compare_string2, $i, 1);
	
	my $add = $c1 + $c2;
	
	$add_string = $add_string.$add;
	
    }
    
    return($add_string, $max, $min);
    
}

sub check_exon_overlap {

    my($gtf_file, $test_start, $test_end) = @_;


    my $gtf = GTF::new({gtf_filename => $gtf_file});
    my $genes = $gtf->transcripts;
    
    foreach my $gene (@$genes) {
	
	my $cds = $gene->cds;
	foreach my $c (@$cds) {
	    


	    if($c->start  < $test_start && $c->end < $test_start) {
	    }
	    elsif($c->start > $test_end && $c->end > $test_end) {
	    }
	    else {
		return("yes");
	    }

	    

	    
	}
	return("no");
    }
    
}

sub build_gtf_filenames {

    my($logfile) = @_;

    # build the names of the gtf files you need to compare
    
    my $thresh;
    if($logfile =~ /.iter_log\.(\d+)$/) {
	$thresh = $1;
    }
    else {
	die("improperly formatted logfile name");
    }
    
    my $gtf1 = $logfile;
    if($gtf1 =~ /\.iter_log\.\d+$/) {
	$gtf1 =~ s/iter_log/first_gtf2/;
    }
    else {
	die("improperly formatted logfile name");
    }
    
    my $gtf2 = $logfile;
    if($gtf2 =~ /\.iter_log\.\d+$/) {
	$gtf2 =~ s/iter_log/second_gtf2/;
    }
    else {
	die("improperly formatted logfile name");
    }
    

    return($thresh, $gtf1, $gtf2);

}


sub do_5_prime_splice_site {
    
    my ($altsplice_bouandries, $one_past_gtf1, $one_past_gtf2, $representative_threshold, $most_numerous_key) = @_;

    my $gtf1_exon1_overlap = check_exon_overlap($one_past_gtf1, $altsplice_bouandries->{$representative_threshold}->{$most_numerous_key}->{0}, $altsplice_bouandries->{$representative_threshold}->{$most_numerous_key}->{2});
    my $gtf2_exon1_overlap = check_exon_overlap($one_past_gtf2, $altsplice_bouandries->{$representative_threshold}->{$most_numerous_key}->{0}, $altsplice_bouandries->{$representative_threshold}->{$most_numerous_key}->{2});
    my $gtf1_exon2_overlap = check_exon_overlap($one_past_gtf1, $altsplice_bouandries->{$representative_threshold}->{$most_numerous_key}->{3}, $altsplice_bouandries->{$representative_threshold}->{$most_numerous_key}->{4});
    my $gtf2_exon2_overlap = check_exon_overlap($one_past_gtf2, $altsplice_bouandries->{$representative_threshold}->{$most_numerous_key}->{3}, $altsplice_bouandries->{$representative_threshold}->{$most_numerous_key}->{4});


    if($gtf1_exon1_overlap eq "yes" &&
       $gtf2_exon1_overlap eq "no" &&
       $gtf1_exon2_overlap eq "yes" &&
       $gtf2_exon2_overlap eq "yes") {

	return("do_5_prime_splice_site, exon 1 dissapears");

    }
    elsif($gtf1_exon1_overlap eq "yes" &&
       $gtf2_exon1_overlap eq "yes" &&
       $gtf1_exon2_overlap eq "yes" &&
       $gtf2_exon2_overlap eq "no") {
	return("do_5_prime_splice_site, exon 2 dissapears");
    }
    elsif($gtf1_exon1_overlap eq "yes" &&
       $gtf2_exon1_overlap eq "no" &&
       $gtf1_exon2_overlap eq "yes" &&
       $gtf2_exon2_overlap eq "no") {
	return("do_5_prime_splice_site, both secondary alignment exons dissapear");
    }
    elsif($gtf1_exon1_overlap eq "yes" &&
       $gtf2_exon1_overlap eq "yes" &&
       $gtf1_exon2_overlap eq "yes" &&
       $gtf2_exon2_overlap eq "yes") {
	return("do_5_prime_splice_site, both secondary alignment exons visible");
    }
    else {
	return("do_5_prime_splice_site unexpected case, e-mail author");
    }

    return("do_5_prime_splice_site,  unexpectedly reached end of function e-mail author");

}

sub do_3_prime_splice_site {
    
    my ($altsplice_bouandries, $one_past_gtf1, $one_past_gtf2, $representative_threshold, $most_numerous_key) = @_;

    my $gtf1_exon1_overlap = check_exon_overlap($one_past_gtf1, $altsplice_bouandries->{$representative_threshold}->{$most_numerous_key}->{0}, $altsplice_bouandries->{$representative_threshold}->{$most_numerous_key}->{1});
    my $gtf2_exon1_overlap = check_exon_overlap($one_past_gtf2, $altsplice_bouandries->{$representative_threshold}->{$most_numerous_key}->{0}, $altsplice_bouandries->{$representative_threshold}->{$most_numerous_key}->{1});
    my $gtf1_exon2_overlap = check_exon_overlap($one_past_gtf1, $altsplice_bouandries->{$representative_threshold}->{$most_numerous_key}->{2}, $altsplice_bouandries->{$representative_threshold}->{$most_numerous_key}->{4});
    my $gtf2_exon2_overlap = check_exon_overlap($one_past_gtf2, $altsplice_bouandries->{$representative_threshold}->{$most_numerous_key}->{2}, $altsplice_bouandries->{$representative_threshold}->{$most_numerous_key}->{4});


    if($gtf1_exon1_overlap eq "yes" &&
       $gtf2_exon1_overlap eq "no" &&
       $gtf1_exon2_overlap eq "yes" &&
       $gtf2_exon2_overlap eq "yes") {

	return("do_3_prime_splice_site, exon 1 dissapears");

    }
    elsif($gtf1_exon1_overlap eq "yes" &&
       $gtf2_exon1_overlap eq "yes" &&
       $gtf1_exon2_overlap eq "yes" &&
       $gtf2_exon2_overlap eq "no") {
	return("do_3_prime_splice_site, exon 2 dissapears");
    }
    elsif($gtf1_exon1_overlap eq "yes" &&
       $gtf2_exon1_overlap eq "no" &&
       $gtf1_exon2_overlap eq "yes" &&
       $gtf2_exon2_overlap eq "no") {
	return("do_3_prime_splice_site, both secondary alignment exons dissapear");
    }
    elsif($gtf1_exon1_overlap eq "yes" &&
       $gtf2_exon1_overlap eq "yes" &&
       $gtf1_exon2_overlap eq "yes" &&
       $gtf2_exon2_overlap eq "yes") {
	return("do_3_prime_splice_site, both secondary alignment exons visible");
    }
    else {
	return("do_3_prime_splice_site unexpected case, e-mail authof");
    }

    return("do_3_prime_splice_site,  unexpectedly reached end of function e-mail author");

}

sub do_retained_intron {
    
    my ($altsplice_bouandries, $one_past_gtf1, $one_past_gtf2, $representative_threshold, $most_numerous_key) = @_;
    

    my $gtf1_exon1_overlap = check_exon_overlap($one_past_gtf1, $altsplice_bouandries->{$representative_threshold}->{$most_numerous_key}->{1}, $altsplice_bouandries->{$representative_threshold}->{$most_numerous_key}->{2});
    my $gtf2_exon1_overlap = check_exon_overlap($one_past_gtf2, $altsplice_bouandries->{$representative_threshold}->{$most_numerous_key}->{1}, $altsplice_bouandries->{$representative_threshold}->{$most_numerous_key}->{2});
    my $gtf1_exon2_overlap = check_exon_overlap($one_past_gtf1, $altsplice_bouandries->{$representative_threshold}->{$most_numerous_key}->{3}, $altsplice_bouandries->{$representative_threshold}->{$most_numerous_key}->{4});
    my $gtf2_exon2_overlap = check_exon_overlap($one_past_gtf2, $altsplice_bouandries->{$representative_threshold}->{$most_numerous_key}->{3}, $altsplice_bouandries->{$representative_threshold}->{$most_numerous_key}->{4});



    if($gtf1_exon1_overlap eq "yes" &&
       $gtf2_exon1_overlap eq "no" &&
       $gtf1_exon2_overlap eq "yes" &&
       $gtf2_exon2_overlap eq "yes") {

	return("do_retained_intron, exon 1 dissapears");

    }
    elsif($gtf1_exon1_overlap eq "yes" &&
       $gtf2_exon1_overlap eq "yes" &&
       $gtf1_exon2_overlap eq "yes" &&
       $gtf2_exon2_overlap eq "no") {
	return("do_retained_intron, exon 2 dissapears");
    }
    elsif($gtf1_exon1_overlap eq "yes" &&
       $gtf2_exon1_overlap eq "no" &&
       $gtf1_exon2_overlap eq "yes" &&
       $gtf2_exon2_overlap eq "no") {
	return("do_retained_intron, both secondary alignment exons dissapear");
    }
    elsif($gtf1_exon1_overlap eq "yes" &&
       $gtf2_exon1_overlap eq "yes" &&
       $gtf1_exon2_overlap eq "yes" &&
       $gtf2_exon2_overlap eq "yes") {
	return("do_retained_intron, both secondary alignment exons visible");
    }
    else {
	return("do_retained_intron unexpected case e-mail author");
    }

    return("do_retained_intron,  unexpectedly reached end of function e-mail author");

}

sub do_clean_alternate_exon {
    
    my ($altsplice_bouandries, $one_past_gtf1, $one_past_gtf2, $representative_threshold, $most_numerous_key) = @_;
    
    my $gtf1_exon1_overlap = check_exon_overlap($one_past_gtf1, $altsplice_bouandries->{$representative_threshold}->{$most_numerous_key}->{0}, $altsplice_bouandries->{$representative_threshold}->{$most_numerous_key}->{1});
    my $gtf2_exon1_overlap = check_exon_overlap($one_past_gtf2, $altsplice_bouandries->{$representative_threshold}->{$most_numerous_key}->{0}, $altsplice_bouandries->{$representative_threshold}->{$most_numerous_key}->{1});
    my $gtf1_exon2_overlap = check_exon_overlap($one_past_gtf1, $altsplice_bouandries->{$representative_threshold}->{$most_numerous_key}->{4}, $altsplice_bouandries->{$representative_threshold}->{$most_numerous_key}->{5});
    my $gtf2_exon2_overlap = check_exon_overlap($one_past_gtf2, $altsplice_bouandries->{$representative_threshold}->{$most_numerous_key}->{4}, $altsplice_bouandries->{$representative_threshold}->{$most_numerous_key}->{5});


    if($gtf1_exon1_overlap eq "yes" &&
       $gtf2_exon1_overlap eq "no" &&
       $gtf1_exon2_overlap eq "yes" &&
       $gtf2_exon2_overlap eq "yes") {

	return("do_clean_alternate_exon,  looks like exon 2 is the first double trace exon past the cassette exon");

    }
    elsif($gtf1_exon1_overlap eq "yes" &&
       $gtf2_exon1_overlap eq "yes" &&
       $gtf1_exon2_overlap eq "yes" &&
       $gtf2_exon2_overlap eq "no") {
	return("do_clean_alternate_exon,   looks like exon 1 is the first double trace exon past the cassette exon");
    }
    elsif($gtf1_exon1_overlap eq "yes" &&
       $gtf2_exon1_overlap eq "no" &&
       $gtf1_exon2_overlap eq "yes" &&
       $gtf2_exon2_overlap eq "no") {
	return("do_clean_alternate_exon, both secondary alignment exons dissapear");
    }
    elsif($gtf1_exon1_overlap eq "yes" &&
       $gtf2_exon1_overlap eq "yes" &&
       $gtf1_exon2_overlap eq "yes" &&
       $gtf2_exon2_overlap eq "yes") {
	return("do_clean_alternate_exon, both secondary alignment exons visible");
    }
    else {
	return("do_clean_alternate_exon unexpected case e-mail author");
    }

    return("do_clean_alternate_exon,  unexpectedly reached end of function e-mail author");

}

sub do_alternate_exon {
    
    my ($altsplice_bouandries, $one_past_gtf1, $one_past_gtf2, $representative_threshold, $most_numerous_key) = @_;
    
    my $gtf1_exon1_overlap = check_exon_overlap($one_past_gtf1, $altsplice_bouandries->{$representative_threshold}->{$most_numerous_key}->{0}, $altsplice_bouandries->{$representative_threshold}->{$most_numerous_key}->{2});
    my $gtf2_exon1_overlap = check_exon_overlap($one_past_gtf2, $altsplice_bouandries->{$representative_threshold}->{$most_numerous_key}->{0}, $altsplice_bouandries->{$representative_threshold}->{$most_numerous_key}->{2});
    my $gtf1_exon2_overlap = check_exon_overlap($one_past_gtf1, $altsplice_bouandries->{$representative_threshold}->{$most_numerous_key}->{5}, $altsplice_bouandries->{$representative_threshold}->{$most_numerous_key}->{7});
    my $gtf2_exon2_overlap = check_exon_overlap($one_past_gtf2, $altsplice_bouandries->{$representative_threshold}->{$most_numerous_key}->{5}, $altsplice_bouandries->{$representative_threshold}->{$most_numerous_key}->{7});


    if($gtf1_exon1_overlap eq "yes" &&
       $gtf2_exon1_overlap eq "no" &&
       $gtf1_exon2_overlap eq "yes" &&
       $gtf2_exon2_overlap eq "yes") {

	return("do_alternate_exon,  looks like exon 2 is the first double trace exon past the cassette exon");

    }
    elsif($gtf1_exon1_overlap eq "yes" &&
       $gtf2_exon1_overlap eq "yes" &&
       $gtf1_exon2_overlap eq "yes" &&
       $gtf2_exon2_overlap eq "no") {
	return("do_alternate_exon,   looks like exon 1 is the first double trace exon past the cassette exon");
    }
    elsif($gtf1_exon1_overlap eq "yes" &&
       $gtf2_exon1_overlap eq "no" &&
       $gtf1_exon2_overlap eq "yes" &&
       $gtf2_exon2_overlap eq "no") {
	return("do_alternate_exon, both secondary alignment exons dissapear");
    }
    elsif($gtf1_exon1_overlap eq "yes" &&
       $gtf2_exon1_overlap eq "yes" &&
       $gtf1_exon2_overlap eq "yes" &&
       $gtf2_exon2_overlap eq "yes") {
	return("do_alternate_exon, both secondary alignment exons visible");
    }
    else {
	return("do_alternate_exon unexpected case, e-mail author");
    }

    return("do_alternate_exon,  unexpectedly reached end of function, e-mail author");

}

sub analyze_all_altsplices {

    my($alt_type_hash, $movie_log, $regex_type) = @_;

    
    # set the regular expression corresponding to the type of alternate
    # splice the function was called to handle
    my $regex = "";
    my $landmark_count = -1;
    if($regex_type eq "alternate exon") {
	$regex = "(2+)(1*)(0+)(1+)(0+)(1*)(2+)";
	$landmark_count = 7;
    }
    elsif($regex_type eq "clean alternate exon") {
	$regex = "(2+)(0+)(1+)(0+)(2+)";
	$landmark_count = 5;
    }
    elsif($regex_type eq "alt 3\' splice site") {
	$regex = "(2+)(0+)(1+)(2+)";
	$landmark_count = 4;
    }
    elsif($regex_type eq "alt 5\' splice site") {
	$regex = "(2+)(1+)(0+)(2+)";
	$landmark_count = 4;
    }
    elsif($regex_type eq "retained intron") {
	$regex = "(0+)(2+)(1+)(2+)(1*)\$|^(1*)(2+)(1+)(2+)(0+)|^(1*)(2+)(1+)(2+)(1*)\$|(0+)(2+)(1+)(2+)(0+)";
	$landmark_count = 5;
    }
    else {
	die("analyze_all altsplices:  error, no test defined for $regex_type");
    }
    
    # open the movie log file
    open(OUT, ">>${movie_log}");

    my %compare_altsplice_hash;
    my %altsplice_bouandries;

    # look at each logfile that matches the type of alternate splice you are looking at
    foreach my $logfile (@{$alt_type_hash->{$regex_type}}) {

        # get the threshold and gtf filenames you will need to use
	my($thresh, $gtf1, $gtf2) = build_gtf_filenames($logfile);

        # added this to fix a problem that occurs when you next out of the loop you need that is one past the last observed alternate splice
	if(!(-e $gtf1) || !(-e $gtf2)) {
	    print OUT "module report:  one of files $gtf1 or $gtf2 needed for analysis doesn't exits... exiting\n";
	    exit(0);
	}

        # generate the add string and coordinates that orient the add string to the genomic sequence
	my ($add_string, $add_string_max, $add_string_min) = make_add_string($gtf1, $gtf2);

        # look for all places in the add string that match the regular
        # expression for the alternate splice... hopefully there is
        # only one but have implemented code to handle the case where
        # there is more than 1
	my $hit_count = 0;
	while($add_string =~ /$regex/g) {
	    $hit_count++;

	    my @landmark_string_array = ($1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17, $18, $19, $20);

            # since there are 4 different patterns in the retained intron regex
            # now need to account for the fact that you might not match the first
            # one, if not you need to shift things down

	    if($regex_type eq "retained intron") {

		if(defined($1)) {
                   # do nothing... the landmark strings are where they need to be
		}
		elsif(defined($6)) {
		    $landmark_string_array[0] = $landmark_string_array[5];
		    $landmark_string_array[1] = $landmark_string_array[6];
		    $landmark_string_array[2] = $landmark_string_array[7];
		    $landmark_string_array[3] = $landmark_string_array[8];
		    $landmark_string_array[4] = $landmark_string_array[9];
		}
		elsif(defined($11)) {
		    $landmark_string_array[0] = $landmark_string_array[10];
		    $landmark_string_array[1] = $landmark_string_array[11];
		    $landmark_string_array[2] = $landmark_string_array[12];
		    $landmark_string_array[3] = $landmark_string_array[13];
		    $landmark_string_array[4] = $landmark_string_array[14];
		}
		elsif(defined($16)) {
		    $landmark_string_array[0] = $landmark_string_array[15];
		    $landmark_string_array[1] = $landmark_string_array[16];
		    $landmark_string_array[2] = $landmark_string_array[17];
		    $landmark_string_array[3] = $landmark_string_array[18];
		    $landmark_string_array[4] = $landmark_string_array[19];
		}
		else {
		    die("error when testing for retained intron, exiting");
		}

	    }

	    my @landmark_length_array = ();
	    for(my $i=0;$i<$landmark_count;$i++) {
		$landmark_length_array[$i] = length($landmark_string_array[$i]);
	    }


	    my $pos = pos($add_string);
	    for(my $i=0;$i<$landmark_count;$i++) {
		$pos -= $landmark_length_array[$i];
	    }

	    my @landmark_array = ();

	    $landmark_array[0] = $add_string_min + $pos;
	    for(my $i=1;$i<$landmark_count+1;$i++) {
		$landmark_array[$i] = $landmark_array[$i-1] + $landmark_length_array[$i-1];
	    }


	    my $compare_altsplice_key;
	    if($regex_type eq "alternate exon") {
		$compare_altsplice_key = $landmark_array[1]."_".$landmark_array[2]."_".$landmark_array[3]."_".$landmark_array[4]."_".$landmark_array[5];
	    }
	    elsif($regex_type eq "clean alternate exon") {
		$compare_altsplice_key = $landmark_array[1]."_".$landmark_array[2]."_".$landmark_array[3]."_".$landmark_array[4];
	    }
	    elsif($regex_type eq "alt 3\' splice site") {
		$compare_altsplice_key = $landmark_array[1]."_".$landmark_array[2]."_".$landmark_array[3];
	    }
	    elsif($regex_type eq "alt 5\' splice site") {
		$compare_altsplice_key = $landmark_array[1]."_".$landmark_array[2]."_".$landmark_array[3];
	    }
	    elsif($regex_type eq "retained intron") {
		$compare_altsplice_key = $landmark_array[1]."_".$landmark_array[2]."_".$landmark_array[3]."_".$landmark_array[4];
	    }
	    else {
		die("attempt to test for unexpected alternate splice type $regex_type... exiting");
	    }


	    push(@{$compare_altsplice_hash{$compare_altsplice_key}}, $logfile);

	    for(my $i=0;$i<$landmark_count+1;$i++) {
		$altsplice_bouandries{$thresh}{$compare_altsplice_key}{$i} = $landmark_array[$i];
	    }
	    
	}

        if($hit_count == 0) {
	    die("should have found alternate splice in $logfile but didnt... exiting, e-mail author");
	}

    }

    my $num_altsplices = keys (%compare_altsplice_hash);

    print OUT "saw $num_altsplices different example(s) of this type of alternate splice\n";
 

    my $most_numerous = -1;
    my $most_numerous_key = "";
    
    if($num_altsplices == 0) {
	return("should have seen at least one altsplice of type $regex_type here but didnt...");
    }

    if($num_altsplices == 1) {
	$most_numerous_key = (keys %compare_altsplice_hash)[0];
	$most_numerous = @{$compare_altsplice_hash{$most_numerous_key}};
    }
    elsif($num_altsplices > 1) {

    # only analyze the most numerous form and check to make sure that it
    # hits on at least 2 thresholds... otherwise they tend to look junky
    # like alignment artifects... will loose some good ones this way but
    # will probably loose a lot more crumy ones

	foreach my $key (keys %compare_altsplice_hash) {

	    my $n = @{$compare_altsplice_hash{$key}};

	    if($n > $most_numerous) {
		$most_numerous = $n;
		$most_numerous_key = $key;
	    }

	    print OUT $key, "\n";
	    foreach my $thing (@{$compare_altsplice_hash{$key}}) {
		print OUT "\t", $thing, "\n";
	    }


	}

        # delete the entries corresponding to the less numerous alternate
        # splice forms seen and proceed as before with the most numerous one 
	my @delete_list = ();
	foreach my $key (keys %compare_altsplice_hash) {
	    if($key ne $most_numerous_key) {
		push(@delete_list, $key)
	    }
	}


	foreach my $key (@delete_list) {
	    delete $compare_altsplice_hash{$key};
	}

	
    }


    my $min_thresh = 1000;
    my $max_thresh = -1;

    foreach my $key (keys %compare_altsplice_hash) {

	foreach my $logfile (@{$compare_altsplice_hash{$key}}) {

	    my $alt_thresh;

	    if($logfile =~ /.iter_log\.(\d+)/) {
		$alt_thresh = $1;
	    }
	    else {
		die("improperly formatted logfile name:\n\n $logfile\n\n");
	    }

	    if($min_thresh > $alt_thresh) {
		$min_thresh = $alt_thresh
	    }

	    if($max_thresh < $alt_thresh) {
		$max_thresh = $alt_thresh;
	    }

	}
	print OUT "looks like range on which I see the alternate splice is $min_thresh to $max_thresh...\n";
    }

# now have the range of thresholds on which we can see the alternate
# splice... go one past the max thresh and make sure one of the
# flanking exons (the double trace side) overlaps in both TR
# alignments but that other flanking exon (the now masked-by-noise
# single trace side) doesnt.


    if($max_thresh == 20) {
	return("$regex_type, see alternate splice all the way to threshold of 20... wont be able to see any double trace only aligments!");
    }


    # get the gtf files for the alignments one past the last alternate
    # splice you saw to check for a double trace only alignment


    my $thresh_to_test = $max_thresh + 1;

    my $one_past_gtf1 = $movie_log;
    $one_past_gtf1 =~ s/movie_log/first_gtf2/;
    $one_past_gtf1 .= ".${thresh_to_test}";
    my $one_past_gtf2 = $one_past_gtf1;
    $one_past_gtf2 =~ s/first_gtf2/second_gtf2/;

    my $temp_logfile_for_threshold = ${$compare_altsplice_hash{$most_numerous_key}}[0];

    my $representative_threshold;
    if($temp_logfile_for_threshold =~ /\.iter_log\.(\d+)/) {
	$representative_threshold = $1;
    }
    else {
	die("error... e-mail author");
    }


    if($regex_type eq "alternate exon") {
	my $return_string = do_alternate_exon(\%altsplice_bouandries, $one_past_gtf1, $one_past_gtf2, $representative_threshold, $most_numerous_key);
	return($return_string);
    }
    elsif($regex_type eq "clean alternate exon") {
	my $return_string = do_clean_alternate_exon(\%altsplice_bouandries, $one_past_gtf1, $one_past_gtf2, $representative_threshold, $most_numerous_key);
	return($return_string);
    }
    elsif($regex_type eq "alt 3\' splice site") {
	my $return_string = do_3_prime_splice_site(\%altsplice_bouandries, $one_past_gtf1, $one_past_gtf2, $representative_threshold, $most_numerous_key);
	return($return_string);
    }
    elsif($regex_type eq "alt 5\' splice site") {
	my $return_string = do_5_prime_splice_site(\%altsplice_bouandries, $one_past_gtf1, $one_past_gtf2, $representative_threshold, $most_numerous_key);
	return($return_string);
    }
    elsif($regex_type eq "retained intron") {
	my $return_string = do_retained_intron(\%altsplice_bouandries, $one_past_gtf1, $one_past_gtf2, $representative_threshold, $most_numerous_key);
	return($return_string);
    }
    else {
	return("analyze_all altsplices:  no test defined for $regex_type");
    }

}


sub analyze_movie {
    my($logfile_list, $movie_log) = @_;

    my %alt_type_hash;
    foreach my $logfile (@{$logfile_list}) {

	open(IN, $logfile);
	chomp(my @in = <IN>);
	close(IN);
	my $result = $in[0];
	
	my $type;
	
	if($result =~ "no altsplice found") {
	    $type = "none";
	    push(@{$alt_type_hash{$type}}, $logfile);
	}
	elsif($result =~ /second alignment was not foward to foward/) {
	    
	}
	else {
	    my @line = split(/\t/, $result);
	    $type = $line[1];

	    if($type =~ /alternate 3\' splice site/) {
		$type = "alt 3' splice site";
	    }
	    elsif($type =~ /alternative 5\' splice site/) {
		$type = "alt 5' splice site";
	    }
	    
	    push(@{$alt_type_hash{$type}}, $logfile);
	}
	
    }
    
    my $num_types_seen = (keys %alt_type_hash);
    open(OUT, ">${movie_log}");
    print OUT "number of types seen = $num_types_seen\n";
    foreach my $thing (keys %alt_type_hash) {
	my $n = @{$alt_type_hash{$thing}};
	print OUT $thing, "\t", $n, "\n";
    }
    
    my $module_report = "null";
    if($num_types_seen == 1) {
	my $only_type = (keys %alt_type_hash)[0];
	if($only_type eq "none") {
	    $module_report = "no alternate splices detected in any of the thresholds\n";
	}
	else {
	    $module_report = "alternate splice of type $only_type seen in all thresholds... error, e-mail author";
	}

    }
    elsif($num_types_seen == 2) {

	print OUT "only one type of alternate splice seen\n";

         # figure out which type youre dealing with and call subroutines to analyze each one

	my $type;
	foreach my $key (keys %alt_type_hash) {
	    if($key ne "none") {
		$type = $key;
	    }
	}


	if($type eq "clean alternate exon") {
	    $module_report = analyze_all_altsplices(\%alt_type_hash, $movie_log, "clean alternate exon");
	}
	elsif($type eq "alt 3\' splice site") {
	    $module_report = analyze_all_altsplices(\%alt_type_hash, $movie_log, "alt 3\' splice site");
	}
	elsif($type eq "alt 5\' splice site") {
	    $module_report = analyze_all_altsplices(\%alt_type_hash, $movie_log, "alt 5\' splice site");
	}
	elsif($type eq "retained intron") {
	    $module_report = analyze_all_altsplices(\%alt_type_hash, $movie_log, "retained intron");
	}
	elsif($type eq "alternate exon") {
	    $module_report = analyze_all_altsplices(\%alt_type_hash, $movie_log, "alternate exon");
	}
	else {
	    $module_report =  "not handling alternate splice type $type, error";
	}

    }
    else {
         # will have to find most common one here and analyze as above

	my $most_common_count = -1;
	my $most_common_type = "not defined";
	foreach my $key (keys %alt_type_hash) {

	    if($key ne "none") {

		my $n = @{$alt_type_hash{$key}};
		if($n > $most_common_count) {
		    $most_common_count = $n;
		    $most_common_type = $key;
		}

	    }

	}

	my %temp_alt_type_hash;

	foreach my $thing (@{$alt_type_hash{"none"}}) {
	    push(@{$temp_alt_type_hash{"none"}}, $thing);
	}

	foreach my $thing (@{$alt_type_hash{$most_common_type}}) {
	    push(@{$temp_alt_type_hash{$most_common_type}}, $thing);
	}

	if($most_common_type eq "clean alternate exon") {
	    $module_report = analyze_all_altsplices(\%alt_type_hash, $movie_log, "clean alternate exon");
	}
	elsif($most_common_type eq "alt 3\' splice site") {
	    $module_report = analyze_all_altsplices(\%alt_type_hash, $movie_log, "alt 3\' splice site");
	}
	elsif($most_common_type eq "alt 5\' splice site") {
	    $module_report = analyze_all_altsplices(\%alt_type_hash, $movie_log, "alt 5\' splice site");
	}
	elsif($most_common_type eq "retained intron") {
	    $module_report = analyze_all_altsplices(\%alt_type_hash, $movie_log, "retained intron");
	}
	elsif($most_common_type eq "alternate exon") {
	    $module_report = analyze_all_altsplices(\%alt_type_hash, $movie_log, "alternate exon");
	}
	else {
	    $module_report = "not handling type $most_common_type, error";
	}

    }
    
    open(OUT, ">>${movie_log}");
    print OUT "module report:  $module_report\n";
    return;
}




sub read_fasta {
    my($filename, $hash_ref) = @_;
    open(FASTA, $filename) || die "couldnt open fasta file $filename\n";
    chomp(my @fasta = <FASTA>);

    my $read_name;
    foreach my $line (@fasta) {
	if($line =~ /^>(.*)/) {
	    $read_name = $1;
	}
	else {
	    $hash_ref->{$read_name} .= $line;
	}
    }
    close(FASTA);   

}

sub return_fasta_seq {

    my ($filename) = @_;

    my %hash;
    read_fasta($filename, \%hash);
    my $key = (keys %hash)[0];
    my $sequence = uc($hash{$key});

    return($sequence);
}

sub populate_map_hash {
    my ($map_file, $map_hash) = @_;

    open(IN, $map_file);
    chomp(my @in = <IN>);
    close(IN);

    foreach my $line (@in) {

	my @line = split(/\s+/, $line);

	if($line[0] ne "refseq1") {

	    $map_hash->{$line[1]} = $line[0];

	}

    }

}

sub get_evidence_list {

    my ($evidence_file) = @_;
    
    open(IN, $evidence_file);
    chomp(my @evidence_list = <IN>);
    close(IN);

    return(@evidence_list);

}


sub guided_generate_sequences {

    my ($poly_hash) = @_;
    
    my $new_ambig = "";
    my $new_recall = "";

    for(my $i=0;$i<$poly_hash->{"length"};$i++) {

	if($poly_hash->{$i}->{"base1_dest"} eq "?" || $poly_hash->{$i}->{"base2_dest"} eq "?") {
	    $new_ambig .= "X";
	    $new_recall .= "X";
	}
	elsif($poly_hash->{$i}->{"base1_dest"} eq "both" || $poly_hash->{$i}->{"base2_dest"} eq "neither") {
	    $new_ambig  .= $poly_hash->{$i}->{"base1"};
	    $new_recall .= $poly_hash->{$i}->{"base1"};
	}
	elsif($poly_hash->{$i}->{"base1_dest"} eq "neither" || $poly_hash->{$i}->{"base2_dest"} eq "both") {
	    $new_ambig  .= $poly_hash->{$i}->{"base2"};
	    $new_recall .= $poly_hash->{$i}->{"base2"};
	}
	elsif($poly_hash->{$i}->{"base1_dest"} eq "ambig" || $poly_hash->{$i}->{"base2_dest"} eq "recall") {
	    $new_ambig  .= $poly_hash->{$i}->{"base1"};
	    $new_recall .= $poly_hash->{$i}->{"base2"};
	}
	elsif($poly_hash->{$i}->{"base1_dest"} eq "recall" || $poly_hash->{$i}->{"base2_dest"} eq "ambig") {
	    $new_ambig  .= $poly_hash->{$i}->{"base2"};
	    $new_recall .= $poly_hash->{$i}->{"base1"};
	}
	elsif($poly_hash->{$i}->{"base1_dest"} eq "ambig" || $poly_hash->{$i}->{"base2_dest"} eq "nowhere") {
	    $new_ambig  .= $poly_hash->{$i}->{"base1"};
	    $new_recall .= "X";
	}
	elsif($poly_hash->{$i}->{"base2_dest"} eq "ambig" || $poly_hash->{$i}->{"base1_dest"} eq "nowhere") {
	    $new_ambig  .= $poly_hash->{$i}->{"base2"};
	    $new_recall .= "X";
	}
	else {
	    die("encountered an unexpected error in function guided_generate_sequences");
	}
	

    }

    return($new_ambig, $new_recall);

}


sub classify_position {

    my($type_hash_ptr, $first_est_char, $first_symbol_char, $second_symbol_char, $first_genome_char, $second_genome_char) = @_;

    my $tag = "unclassified";

# this is a check to make sure you have both alignments to work with
    if($first_genome_char ne "" && $first_genome_char ne "unaligned" && $second_genome_char ne "" && $second_genome_char ne "unaligned") {

	my $ambig_char_type = $type_hash_ptr->{$first_est_char};

	if($ambig_char_type eq "single") {
	    
	    if($first_symbol_char eq "|") {
		
		if($second_symbol_char eq "|") {
		    # get this tag if you are a single peak and you match in both alignments
		    $tag = "correct single peak";
		}
		elsif($second_symbol_char eq " ") {
		    # get this tag if you are a single peak and you match in the first but dont match in the second
		    # (ie there was a secondary peak that didnt get called... not actually checking for that though)
		    $tag = "missed double peak";
		}
		else {
		    die("not handling this second symbol char type $second_symbol_char");
		}
		
	    }
	    elsif($first_symbol_char eq " ") {
		$tag = "totally blew position";
	    }
	    else {
		die("not handling this first symbol char type $first_symbol_char");
	    }
	    
	}
	elsif($ambig_char_type eq "double") {
	    
	    if($first_symbol_char eq "|") {
		
		if($second_symbol_char eq "|") {
		    # get this tag if you are a double peak that matches in both alignments
		    $tag = "correct double peak";
		}
		elsif($second_symbol_char eq " ") {
		    
		    if(uc($first_genome_char) eq uc($second_genome_char)) {
			# get this tag if you are a double peak that
			# matches in the first alignment but
			# mismatches in the second alignment and if
			# you HADN'T been called as a double peak you
			# would have been a match (ie... the first
			# genome character, which matched, is the same
			# as the second)
			$tag = "noise peak";
		    }
		    else {
			# get this tag if you are a double peak that
			# matched in the first, didnt match in the
			# second and wouldnt have matched if you were
			# called as a single peak (there is a
			# secondary peak but it wouldnt have mattered
			# if it got called becasue it was the wrong
			# one)
			$tag = "incorrect double peak";
		    }
		    
		}
		else {
		    die("not handling this second symbol char type $second_symbol_char");
		}
		
	    }
	    elsif($first_symbol_char eq " ") {
		$tag = "totally blew position";
	    }
	    else {
		die("not handling this first symbol char type $first_symbol_char");
	    }
	    
	}
	elsif($ambig_char_type eq "N") {
	    $tag = "totally blew position";
	}
	else {
	    die("not handling this ambig char type $ambig_char_type");
	}
	
	
    }
    # this says that we are in a part of the trace that had a first alignment but no second alignment
    elsif($first_genome_char ne "" && $first_genome_char ne "unaligned" && ($second_genome_char eq "" || $second_genome_char eq "unaligned")) {

	if($first_symbol_char eq "|") {
	    $tag = "correct primary only";
	}
	elsif($first_symbol_char eq " ") {
	    $tag = "totally blew position";
	}
	else {
	    die("not handling thay type of first symbol char $first_symbol_char");
	}
    }
    
    return($tag);
    
}

sub gather_statistics {

    my ($array_ptr) = @_;

    my $sum = 0;
    my $squared_sum = 0;
    my $count = @$array_ptr;

    my $max = -100;
    my $min = 100;

    foreach my $ratio (@$array_ptr) {
	$sum += $ratio;
	$squared_sum += $ratio ** 2;

	if($ratio < $min) {
	    $min = $ratio;
	}

	if($ratio > $max) {
	    $max = $ratio;
	}

    }

    my $mean;
    my $std_dev;

    if($count != 0) {
	
	$mean = $sum / $count;
	$std_dev = sqrt(($squared_sum / $count) - $mean ** 2);
    }
    else {
	$mean = "NA";
	$std_dev = "NA";
    }
    
    
    my $median = $array_ptr->[int($count / 2)];

    return($count, $mean, $std_dev, $median, $max, $min);

}


# little subroutine to check how many known noise and correct double peaks will get misclassified with any given threshold
sub check_for_misclassification {

    my ($threshold_guess, $noise_peaks, $correct_double_peaks) = @_;

    
    my $bad_noise_peaks = 0;
    my $bad_double_peaks = 0;

    foreach my $ratio (@$noise_peaks) {
	if($ratio < $threshold_guess) {
	    $bad_noise_peaks++;
	    print "would miscall this ratio $ratio as a double peak with threshold $threshold_guess\n";
	}
    }
    
    foreach my $ratio (@$correct_double_peaks) {
	if($ratio > $threshold_guess) {
	    $bad_double_peaks++;
	    print "would miscall this ratio $ratio as a noise peak with threshold $threshold_guess\n";
	}
    }

    my $bad_peaks = $bad_noise_peaks + $bad_double_peaks;

    print "with this threshold:\n";
    print "$bad_noise_peaks bad noise peaks\n";
    print "$bad_double_peaks bad double peaks\n";
    print "$bad_peaks total bad peaks\n";


}

# this is a horrible name for this subroutine... it does a bunch of
# things including trying to classify each of the positions in the
# poly file and telling which sequence each of the called bases should
# go into... ambig or recalled, on the next iteration
sub tester {

    my ($poly_hash, $trace_to_first_align, $first_align_contents, $recall_to_second_align, $second_align_contents, $type_hash, $thresh) = @_;

    my @noise_peaks;
    my @correct_double_peaks;

    for(my $i=0;$i<$poly_hash->{"length"};$i++) {

	my $first_align_pos = $trace_to_first_align->{$i};
	
	my $first_genome_char = "";
	my $first_symbol_char = "";
	my $first_est_char    = "";

	if($first_align_pos ne "unaligned") {
	    
	    $first_genome_char = $first_align_contents->{$first_align_pos}->{"genome"};
	    $first_symbol_char = $first_align_contents->{$first_align_pos}->{"symbol"};
	    $first_est_char    = $first_align_contents->{$first_align_pos}->{"est"};
	}
	
	my $second_align_pos = $recall_to_second_align->{$i};

	my $second_genome_char = "";
	my $second_symbol_char = "";
	my $second_est_char    = "";

	if($second_align_pos ne "unaligned") {
	    $second_genome_char = $second_align_contents->{$second_align_pos}->{"genome"};
	    $second_symbol_char = $second_align_contents->{$second_align_pos}->{"symbol"};
	    $second_est_char    = $second_align_contents->{$second_align_pos}->{"est"};
	}

       # try to classify each position of the poly file as a good
       # double peak, good single peak, single peak plus noise, bad
       # double peak etc...
	my $tag = "unassigned";
	$tag = classify_position(\%type_hash, $first_est_char, $first_symbol_char, $second_symbol_char, $first_genome_char, $second_genome_char);
	$poly_hash->{$i}->{"tag"} = $tag;

        # originally I used the base1_dest and base2_dest variables here to
        # sort out the peaks into ambig or recalled sequence but try doing
        # this by setting the position specific threshold of this
        # position... this is kind of doing the same thing by side effect but
        # it is actually what I had originally envisioned this to be... if you
        # set the threshold a little lower than the actual ratio for the
        # position

	if($tag eq "noise peak") {
	    push(@noise_peaks, $poly_hash->{$i}->{"area_ratio"});

            # figure out which of the 2 peaks ISNT the noise peak and never look at the other (noise) peak again
	    if(uc($poly_hash->{$i}->{"base1"}) eq uc($first_genome_char) && uc($poly_hash->{$i}->{"base1"}) eq uc($second_genome_char)) {
		$poly_hash->{$i}->{"base1_dest"} = "both";
		$poly_hash->{$i}->{"base2_dest"} = "neither";
	    }
	    elsif(uc($poly_hash->{$i}->{"base2"}) eq uc($first_genome_char) && uc($poly_hash->{$i}->{"base2"}) eq uc($second_genome_char)) {
		$poly_hash->{$i}->{"base1_dest"} = "neither";
		$poly_hash->{$i}->{"base2_dest"} = "both";
	    }
	    else {
		die("unexpected error... exiting");
	    }

	}
	elsif($tag eq "correct double peak") {
	    push(@correct_double_peaks, $poly_hash->{$i}->{"area_ratio"});

            # figure out which base matched in which alignment and lock them in
	    if(uc($poly_hash->{$i}->{"base1"}) eq uc($first_genome_char) && uc($poly_hash->{$i}->{"base2"}) eq uc($second_genome_char)) {
		$poly_hash->{$i}->{"base1_dest"} = "ambig";
		$poly_hash->{$i}->{"base2_dest"} = "recall";
	    }
	    elsif(uc($poly_hash->{$i}->{"base1"}) eq uc($second_genome_char) && uc($poly_hash->{$i}->{"base2"}) eq uc($first_genome_char)) {
		$poly_hash->{$i}->{"base1_dest"} = "recall";
		$poly_hash->{$i}->{"base2_dest"} = "ambig";
	    }
	    else {
		die("unexpected error... exiting");
	    }

	}
	elsif($tag eq "incorrect double peak") {
	    $poly_hash->{$i}->{"base1_dest"} = "ambig";
	    $poly_hash->{$i}->{"base2_dest"} = "recall";
	}
	elsif($tag eq "totally blew position") {

	}
	elsif($tag eq "missed double peak") {
	    $poly_hash->{$i}->{"base1_dest"} = "both";
	    $poly_hash->{$i}->{"base2_dest"} = "neither";
	}
	elsif($tag eq "correct single peak") {

	    if(uc($poly_hash->{$i}->{"base1"}) eq uc($first_genome_char) && uc($poly_hash->{$i}->{"base1"}) eq uc($second_genome_char)) {
		$poly_hash->{$i}->{"base1_dest"} = "both";
		$poly_hash->{$i}->{"base2_dest"} = "neither";
	    }
	    else {
		die("unexpected error... exiting");
	    }

	}
	elsif($tag eq "unclassified") {

	}
	elsif($tag eq "unassigned") {

	}
	elsif($tag eq "correct primary only") {

	    if(uc($poly_hash->{$i}->{"base1"}) eq uc($first_genome_char)) {
		$poly_hash->{$i}->{"base1_dest"} = "ambig";
		$poly_hash->{$i}->{"base2_dest"} = "nowhere";
	    }
	    elsif(uc($poly_hash->{$i}->{"base2"}) eq uc($first_genome_char)) {
		$poly_hash->{$i}->{"base1_dest"} = "nowhere";
		$poly_hash->{$i}->{"base2_dest"} = "ambig";
	    }
	    else {
		die("unexpected error... exiting");
	    }

	}
	else {
	    die("unexpected tag enountered $tag exiting");
	}

    }

}



sub create_poly_hash {
    
    my ($poly, $poly_hash, $threshold) = @_;

    open(IN, $poly);
    chomp(my @poly_file = <IN>);
    close(IN);
    
    
    
    my $idx = 0;
    foreach my $line (@poly_file) {
	if($idx != 0) {
	    my @line = split(/\s+/, $line);
	    
	    $poly_hash->{$idx-1}->{"base1"} = $line[0];
	    $poly_hash->{$idx-1}->{"area1"} = $line[2];
	    $poly_hash->{$idx-1}->{"base2"} = $line[4];
	    $poly_hash->{$idx-1}->{"area2"} = $line[6];
	    $poly_hash->{$idx-1}->{"base1_dest"} = "?";
	    $poly_hash->{$idx-1}->{"base2_dest"} = "?";
	    $poly_hash->{$idx-1}->{"tag"} = "?";
	    $poly_hash->{$idx-1}->{"thresh"} = $threshold;


	    if($poly_hash->{$idx-1}->{"base2"} eq "N") {
		$poly_hash->{$idx-1}->{"area_ratio"} = "NA";
	    }
	    else {
		my $larger;
		my $smaller;

		if($poly_hash->{$idx-1}->{"area1"} > $poly_hash->{$idx-1}->{"area2"}) {
		    $larger = $poly_hash->{$idx-1}->{"area1"};
		    $smaller = $poly_hash->{$idx-1}->{"area2"};
		}
		else {
		    $larger = $poly_hash->{$idx-1}->{"area2"};
		    $smaller = $poly_hash->{$idx-1}->{"area1"};
		}

		$poly_hash->{$idx-1}->{"area_ratio"} = $larger / $smaller;
	    }
	    
            # this starts out as a '?' and if an alignment locks the position it will change
	    $poly_hash->{$idx-1}->{"lock"} = "?";  
	    

	}
	$idx++;
    }
    
    $poly_hash->{"length"} = $idx - 1;

}

sub dump_poly_hash {

    my($poly_hash) = @_;

    for(my $i=0;$i<$poly_hash->{"length"};$i++) {

	print $i, "\t", $poly_hash->{$i}->{"base1"}, "\t", $poly_hash->{$i}->{"base2"}, "\n";

    }

}

sub blindly_generate_ambig_sequence {
    my($poly_hash, $ambig_hash) = @_;
    my $ambig_seq = "";

    for(my $i=0;$i<$poly_hash->{"length"};$i++) {

        # check ratio threshold for this position... if ratio < threshold call ambiguity base...
	if($poly_hash->{$i}->{"area_ratio"} ne "NA" && 
	   $poly_hash->{$i}->{"area_ratio"}  < $poly_hash->{$i}->{"thresh"}) {
	    
	    if(defined( my $ambig_char = $ambig_hash->{$poly_hash->{$i}->{"base1"}}{$poly_hash->{$i}->{"base2"}})) {
		$ambig_seq .= $ambig_char;
	    }
	    else {
		die("undefined pair in poly hash, exiting e-mail author");
	    }
	}
        # ...otherwise just use the first base
	else {
	    $ambig_seq .= $poly_hash->{$i}->{"base1"};
	}
	
	

    }
    return $ambig_seq;
}


sub populate_type_hash {

    my($type_hash) = @_;

    $type_hash->{"A"} = "single";
    $type_hash->{"a"} = "single";
    $type_hash->{"C"} = "single";
    $type_hash->{"c"} = "single";
    $type_hash->{"G"} = "single";
    $type_hash->{"g"} = "single";
    $type_hash->{"T"} = "single";
    $type_hash->{"t"} = "single";

    $type_hash->{"M"} = "double";
    $type_hash->{"m"} = "double";
    $type_hash->{"R"} = "double";
    $type_hash->{"r"} = "double";
    $type_hash->{"W"} = "double";
    $type_hash->{"w"} = "double";
    $type_hash->{"S"} = "double";
    $type_hash->{"s"} = "double";
    $type_hash->{"K"} = "double";
    $type_hash->{"k"} = "double";
    $type_hash->{"Y"} = "double";
    $type_hash->{"y"} = "double";

    $type_hash->{"N"} = "N";
    $type_hash->{"n"} = "N";

    $type_hash->{"-"} = "gap";
}


sub populate_ambig_hash {

    my($ambig_hash) = @_;
    
    $ambig_hash->{"A"}->{"C"} = "M";
    $ambig_hash->{"C"}->{"A"} = "M";
    
    $ambig_hash->{"A"}->{"G"} = "R";
    $ambig_hash->{"G"}->{"A"} = "R";
    
    $ambig_hash->{"A"}->{"T"} = "W";
    $ambig_hash->{"T"}->{"A"} = "W";
    
    $ambig_hash->{"C"}->{"G"} = "S";
    $ambig_hash->{"G"}->{"C"} = "S";
    
    $ambig_hash->{"C"}->{"T"} = "Y";
    $ambig_hash->{"T"}->{"C"} = "Y";
    
    $ambig_hash->{"G"}->{"T"} = "K";
    $ambig_hash->{"T"}->{"G"} = "K";
    
    $ambig_hash->{"A"}->{"N"} = "A";
    $ambig_hash->{"N"}->{"A"} = "A";
    
    $ambig_hash->{"C"}->{"N"} = "C";
    $ambig_hash->{"N"}->{"C"} = "C";
    
    $ambig_hash->{"G"}->{"N"} = "G";
    $ambig_hash->{"N"}->{"G"} = "G";
    
    $ambig_hash->{"T"}->{"N"} = "T";
    $ambig_hash->{"N"}->{"T"} = "T";
    
    
    $ambig_hash->{"N"}->{"N"} = "N";
    $ambig_hash->{"A"}->{"A"} = "A";
    $ambig_hash->{"C"}->{"C"} = "C";
    $ambig_hash->{"G"}->{"G"} = "G";
    $ambig_hash->{"T"}->{"T"} = "T";
    
}

sub create_alignment_hash {

    my($alignment_file, $alignment_contents) = @_;

    open(INF, $alignment_file);
    my $report;
    if($report = new EGlite(\*INF)){
	# it worked, dont exit
    }
    else {
	die("couldnt open ${alignment_file}!");
    }
    close(INF);

    my($genome, $symbol, $est) = split(/\n/, $report->alignment);

    $alignment_contents->{"est_begin"} = $report->eb  - 1;

    my $genome_cut;
    my $symbol_cut;
    my $est_cut;

# remove any colunm of the alignment that represents an intron
    for(my $i=0;$i<length($genome);$i++) {
	
	my $genome_char = substr($genome, $i, 1);
	my $symbol_char = substr($symbol, $i, 1);
	my $est_char = substr($est, $i, 1);
	
	if($est_char ne ".") {
	    $genome_cut .= $genome_char;
	    $symbol_cut .= $symbol_char;
	    $est_cut .= $est_char;
	}
    }



# create a hash that has all the information in the alignment indexed
# by the colunm in the alignment

    $alignment_contents->{"length"} = length($genome_cut);

    for(my $i=0;$i<$alignment_contents->{"length"};$i++) {
    
	my $genome_char = substr($genome_cut, $i, 1);
	my $symbol_char = substr($symbol_cut, $i, 1);
	my $est_char = substr($est_cut, $i, 1);
	
	$alignment_contents->{$i}->{"genome"} = $genome_char;
	$alignment_contents->{$i}->{"symbol"} = $symbol_char;
	$alignment_contents->{$i}->{"est"} = $est_char;
	
    }

}

sub dump_alignment_hash {

    print "dump_alignment_hash called\n";

    my($alignment_contents) = @_;

    for(my $i=0;$i<$alignment_contents->{"length"};$i++) {

	print $i, "\t", $alignment_contents->{$i}->{"genome"}, " ", $alignment_contents->{$i}->{"symbol"}, " ", $alignment_contents->{$i}->{"est"}, "\n";

    }

    print "alignment length (after removing introns) = ", $alignment_contents->{"length"}, "\n";
    print "est begins aligning at position ", $alignment_contents->{"est_begin"}, "\n";

}

sub create_sequence_alignment_map {

    my($sequence_to_alignment, $alignment_to_sequence, $alignment_contents, $sequence_length) = @_;

    # mark unaligned sequence as unaligned
    for(my $seq_pos=0;$seq_pos<$alignment_contents->{"est_begin"};$seq_pos++) {
	$sequence_to_alignment->{$seq_pos} = "unaligned";
    }

    my $seq_pos = $alignment_contents->{"est_begin"};

    for(my $align_pos=0;$align_pos<$alignment_contents->{"length"};$align_pos++) {
	
	my $est_char = $alignment_contents->{$align_pos}->{"est"};
	my $genome_char = $alignment_contents->{$align_pos}->{"genome"};
	my $symbol_char = $alignment_contents->{$align_pos}->{"symbol"};
	
        # look at alignment colunm and decide what to put in mappings
	if($est_char eq "-") {
	    
	    $alignment_to_sequence->{$align_pos} = "x";
	    
	}
	else {
	    $alignment_to_sequence->{$align_pos} = $seq_pos;
	    $sequence_to_alignment->{$seq_pos} = $align_pos;
	    $seq_pos++;
	}
    }

    for(my $i = $seq_pos;$i<$sequence_length;$i++) {
	$sequence_to_alignment->{$i} = "unaligned";
    }

}

sub dump_sequence_to_alignment_hash {

    print "dump_sequence_to_alignment_hash called\n";

    my($sequence_to_alignment, $poly_file, $alignment_contents, $ambig_hash, $mode) = @_;

    print $poly_file->{"length"}, "\n";

    for(my $seq_pos=0;$seq_pos<$poly_file->{"length"};$seq_pos++) {

	my $align_pos = $sequence_to_alignment->{$seq_pos};

	my $base1 = $poly_file->{$seq_pos}->{"base1"};
	my $base2 = $poly_file->{$seq_pos}->{"base2"};

	my ($genome, $symbol, $est);

	if($align_pos ne "unaligned") {   
	    $genome = $alignment_contents->{$align_pos}->{"genome"};
	    $symbol = $alignment_contents->{$align_pos}->{"symbol"};
	    $est = $alignment_contents->{$align_pos}->{"est"};

            # do some checks on the maps as you dump them out... depends on which map you are constructing 
	    if($mode eq "primary") {
		# if there is an est base there make sure its the right one according to the poly file/ambig sequence
		if($est ne "-") {
		    if($ambig_hash->{$base1}->{$base2} ne $est) {
			die("cant recreate ambig pos");
		    }
		}
	    }
	    elsif($mode eq "secondary") {
		# if there is an est base that actually came from the
		# trace (ie not an N) make sure it matches one of the
		# bases at the poly file at this position
		if($est ne "-" && uc($est) ne "N") {
		    if(uc($est) ne uc($base1) && uc($est) ne uc($base2)) {
			print "est = $est\n";
			print "base1 = $base1\n";
			print "base2 = $base2\n";
			die("recalled base not one of the bases from the poly hash");
		    }
		}
	    }
	    else {
		die("unknown mode");
	    }
	    
	    print $seq_pos, "\t", $align_pos, "\t", $base1, " ", $base2, "\t", $genome, " ", $symbol, " ", $est, "\n";
	}
	else {
	    print $seq_pos, "\t", $align_pos, "\t", $base1, " ", $base2, "\tunaligned\n";
	}




    }

}

sub dump_alignment_to_sequence_hash {

    my($alignment_to_sequence, $alignment_contents) = @_;

    for(my $align_pos=0;$align_pos<$alignment_contents->{"length"};$align_pos++) {

	print $align_pos, "\t", $alignment_to_sequence->{$align_pos}, "\n";

    }

}


sub blindly_generate_recall_sequence {
    my ($poly_hash, $alignment_contents, $sequence_to_alignment, $disambig_hash) = @_;

    my $recall_seq = "";

    # generating the recalled sequence in this way should ensure that it
    # has the same length as the ambiguity sequence and poly file

    for(my $seq_pos=0;$seq_pos<$poly_hash->{"length"};$seq_pos++) {

	my $align_pos = $sequence_to_alignment->{$seq_pos};

	if($align_pos eq "unaligned") {
	    $recall_seq .= "X";
	}
	else {

	    my $genome = $alignment_contents->{$align_pos}->{"genome"};
	    my $symbol = $alignment_contents->{$align_pos}->{"symbol"};
	    my $est    = $alignment_contents->{$align_pos}->{"est"};

	    if(!defined($disambig_hash->{uc($est)}->{uc($genome)})) {
		die("Didnt define what happens when you see a ", uc($est), " and a ", uc($genome), "e-mail author \n");
	    }

	    my $recall_base = $disambig_hash->{uc($est)}->{uc($genome)};

	    $recall_seq .= $recall_base;

	}

    }

    return $recall_seq;

}

sub populate_disambig_hash {
    
    my($disambig_hash) = @_;

    $disambig_hash->{"M"}->{"A"} = "c";
    $disambig_hash->{"M"}->{"C"} = "a";
    $disambig_hash->{"M"}->{"G"} = "n";
    $disambig_hash->{"M"}->{"T"} = "n";
    $disambig_hash->{"M"}->{"N"} = "n";
    
    $disambig_hash->{"R"}->{"A"} = "g";
    $disambig_hash->{"R"}->{"C"} = "n";
    $disambig_hash->{"R"}->{"G"} = "a";
    $disambig_hash->{"R"}->{"T"} = "n";
    $disambig_hash->{"R"}->{"N"} = "n";
    
    $disambig_hash->{"W"}->{"A"} = "t";
    $disambig_hash->{"W"}->{"C"} = "n";
    $disambig_hash->{"W"}->{"G"} = "n";
    $disambig_hash->{"W"}->{"T"} = "a";
    $disambig_hash->{"W"}->{"N"} = "n";
    
    $disambig_hash->{"S"}->{"A"} = "n";
    $disambig_hash->{"S"}->{"C"} = "g";
    $disambig_hash->{"S"}->{"G"} = "c";
    $disambig_hash->{"S"}->{"T"} = "n";
    $disambig_hash->{"S"}->{"N"} = "n";
    
    $disambig_hash->{"Y"}->{"A"} = "n";
    $disambig_hash->{"Y"}->{"C"} = "t";
    $disambig_hash->{"Y"}->{"G"} = "n";
    $disambig_hash->{"Y"}->{"T"} = "c";
    $disambig_hash->{"Y"}->{"N"} = "n";
    
    $disambig_hash->{"K"}->{"A"} = "n";
    $disambig_hash->{"K"}->{"C"} = "n";
    $disambig_hash->{"K"}->{"G"} = "t";
    $disambig_hash->{"K"}->{"T"} = "g";
    $disambig_hash->{"K"}->{"N"} = "n";
    
    
    $disambig_hash->{"A"}->{"A"} = "A";
    $disambig_hash->{"A"}->{"C"} = "A";
    $disambig_hash->{"A"}->{"G"} = "A";
    $disambig_hash->{"A"}->{"T"} = "A";
    $disambig_hash->{"A"}->{"N"} = "A";
    
    $disambig_hash->{"C"}->{"A"} = "C";
    $disambig_hash->{"C"}->{"C"} = "C";
    $disambig_hash->{"C"}->{"G"} = "C";
    $disambig_hash->{"C"}->{"T"} = "C";
    $disambig_hash->{"C"}->{"N"} = "C";
    
    $disambig_hash->{"G"}->{"A"} = "G";
    $disambig_hash->{"G"}->{"C"} = "G";
    $disambig_hash->{"G"}->{"G"} = "G";
    $disambig_hash->{"G"}->{"T"} = "G";
    $disambig_hash->{"G"}->{"N"} = "G";
    
    $disambig_hash->{"T"}->{"A"} = "T";
    $disambig_hash->{"T"}->{"C"} = "T";
    $disambig_hash->{"T"}->{"G"} = "T";
    $disambig_hash->{"T"}->{"T"} = "T";
    $disambig_hash->{"T"}->{"N"} = "T";

    $disambig_hash->{"N"}->{"A"} = "N";
    $disambig_hash->{"N"}->{"C"} = "N";
    $disambig_hash->{"N"}->{"G"} = "N";
    $disambig_hash->{"N"}->{"T"} = "N";
    $disambig_hash->{"N"}->{"N"} = "N";
    
    $disambig_hash->{"-"}->{"A"} = "N";
    $disambig_hash->{"-"}->{"C"} = "N";
    $disambig_hash->{"-"}->{"G"} = "N";
    $disambig_hash->{"-"}->{"T"} = "N";
    
    $disambig_hash->{"A"}->{"-"} = "N";
    $disambig_hash->{"C"}->{"-"} = "N";
    $disambig_hash->{"G"}->{"-"} = "N";
    $disambig_hash->{"T"}->{"-"} = "N";
    $disambig_hash->{"M"}->{"-"} = "N";
    $disambig_hash->{"R"}->{"-"} = "N";
    $disambig_hash->{"W"}->{"-"} = "N";
    $disambig_hash->{"S"}->{"-"} = "N";
    $disambig_hash->{"Y"}->{"-"} = "N";
    $disambig_hash->{"K"}->{"-"} = "N";
    $disambig_hash->{"N"}->{"-"} = "N";
    
}


sub get_estgenome_report {

    my($filename) = @_;

    open(INF, $filename);
    my $report;
    if($report = new EGlite(\*INF)){
    }
    else {
	die("couldnt open $filename");
    }
    
    close(INF);

    return($report);

}
