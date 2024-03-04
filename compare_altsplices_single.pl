#! /usr/bin/perl -w
use strict;
use GTF;
use Getopt::Long;

#my ($ex_file, $compare_file) = @ARGV;
my($gtf_file1, $gtf_file2) = @ARGV;

my($no_both_altsplice) = (0);

GetOptions( "no_both_altsplice"  => \$no_both_altsplice);

my $report_string = "";

my @altsplices_seen;

#open(IN, $compare_file);
#chomp(my @compare_file = <IN>);
#close(IN);
#my %compare_hash;
#foreach my $line (@compare_file) {
#    my ($file1, $file2) = split(/\s+/, $line);
#    $compare_hash{$file1} = $file2;
#}

#open(IN, $ex_file);
#chomp(my @ex_file = <IN>);
#close(IN);

#my $exp = 0;
#my $not_exp = 0;

#my %exp;
#my %not_exp;

#my %hit_hash;

#$gtf_file1 =~ /(\d\d)A\.scf/;

my $trace_number = $1;

#    print $trace_number, "\n";

my $space_flag = 0;

#my $gtf_file2 = $compare_hash{$gtf_file1};

#printf STDERR "comparing $gtf_file1 & $gtf_file2\n";

my $gtf1 = GTF::new({gtf_filename => $gtf_file1});
my $genes1 = $gtf1->transcripts;

my $gtf2 = GTF::new({gtf_filename => $gtf_file2});
my $genes2 = $gtf2->transcripts;

my $max = 0;
my $min = 1000000000000000;

my $length1 = 0;
my $length2 = 0;

$gtf_file1 =~ s/\.scf\.first_gtf//;

foreach my $gene (@$genes1){
    
    my $cds = $gene->cds;
    foreach my $c (@$cds){
#	print $c->start, "\t", $c->stop, "\n";
	
	$length1 += $c->stop - $c->start;
	
	if($c->stop > $max) {
	    $max = $c->stop;
	}
	if($c->start < $min) {
	    $min = $c->start;
	}
	
    }
    
}

#print "\n\n";

foreach my $gene (@$genes2){
    
    my $cds = $gene->cds;
    foreach my $c (@$cds){
#	print $c->start, "\t", $c->stop, "\n";
	$length2 += $c->stop - $c->start;
	
	
	if($c->stop > $max) {
	    $max = $c->stop;
	}
	if($c->start < $min) {
	    $min = $c->start;
	}
    }
}

#print "\n\n";
#print $min, "\t", $max, "\n";
#print "\n\n";

my $compare_string1 = "";
my $compare_string2 = "";
my $add_string = "";

#print $max, "\t", $min, "\n";

for(my $i=0;$i<$max - $min;$i++) {
    
    $compare_string1 .= "0";
    $compare_string2 .= "0";
#    $add_string .= "0";
}

foreach my $gene (@$genes1){
    
    my $cds = $gene->cds;
    foreach my $c (@$cds){
#	print $c->start - $min, "\t", $c->stop - $min, "\n";
	
	my $temp = "";
	for(my $i=0;$i<$c->stop - $c->start;$i++) {
	    $temp .= "1";
	}
	substr($compare_string1, $c->start - $min, $c->stop - $c->start, $temp);
    }
    
}

#print "\n\n";

foreach my $gene (@$genes2){
    
    my $cds = $gene->cds;
    foreach my $c (@$cds){
#	print $c->start - $min, "\t", $c->stop - $min, "\n";
	
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

if($add_string =~ /2+(1*)0+1+0+(1*)2+/) {
    
    my $five_prime_slop = length($1);
    my $three_prime_slop = length($2);
    
#    my $report_string;
    
    if($five_prime_slop == 0 && $three_prime_slop == 0) {
	$space_flag = 1;
	
#	print STDERR $gtf_file1, "\t", "clean alternate exon\t\t";
	$report_string = $gtf_file1."\t"."clean alternate exon\t\t";

	push(@altsplices_seen, "clean alternate exon");

    }
    else {
	$space_flag = 1;
	
#	print STDERR $gtf_file1, "\t", "alternate exon\t\t\t";
	$report_string = $gtf_file1."\t"."alternate_exon\t\t\t";
	push(@altsplices_seen, "alternate exon");
    }
    

    
}

if($add_string =~ /2+(1*)0+1+0+1+0+(1*)2+/) {
    
    my $five_prime_slop = length($1);
    my $three_prime_slop = length($2);
    
#    my $report_string;
    
    if($five_prime_slop == 0 && $three_prime_slop == 0) {
	
	$space_flag = 1;
	
#	print STDERR $gtf_file1, "\t", "clean mutually exclusive exons\t";
	$report_string = $gtf_file1."\t"."clean mutually exclusive exons\t";

	push(@altsplices_seen, "clean mutually exclusive exons");
	    
	
    }
    else {
	$space_flag = 1;
	
#	print STDERR $gtf_file1, "\t", "mutually exclusive exons\t";
	$report_string = $gtf_file1."\t"."mutually exclusive exons\t";	    
	
	push(@altsplices_seen, "mutually exclusive exons");

    }
    
}

if($add_string =~ /2+(1+)0+2+/) {
    
    my $l = length($1);
    
    if($l > 10) {
	
#	my $report_string;
	
	$space_flag = 1;
#	print STDERR $gtf_file1, "\t", "alternate 5' splice site $l\t";
	$report_string = $gtf_file1."\t"."alternate 5' splice site $l\t";

	push(@altsplices_seen, "alternative 5' splice site $l");

	my $flag = 0;
	
    }
    
}

if($add_string =~ /2+0+(1+)2+/) {
    
    my $l = length($1);
    
    if($l > 10) {
	
	$space_flag = 1;
#	print STDERR $gtf_file1, "\t", "alternate 3' splice site $l\t";
	$report_string = $gtf_file1."\t"."alternate 3' splice site $l\t";    
	push(@altsplices_seen, "alternate 3' splice site $l");	
    }
    
}

if(!($no_both_altsplice) && $add_string =~ /2+(1+)0+(1+)2+/) {
    
    my $l1 = length($1);
    my $l2 = length($2);
    
    my $difference = abs($l1 - $l2);

    if(($l1 > 10 || $l2 > 10) && $difference > 10) {
	
	$space_flag = 1;
#	print STDERR $gtf_file1, "\t", "alternate splice sites (both) $l1 $l2\t";
	$report_string = $gtf_file1."\t"."alternate splice sites (both) $l1 $l2\t";	    

	push(@altsplices_seen, "alternate splice sites (both) $l1 $l2");
	
	
    }
    
}


#if($add_string =~ /2+(1+)2+/) {
#    if(length($1) > 10) {
#	$space_flag = 1;
#	$report_string = $gtf_file1."\t"."retained intron\t\t\t";
#	push(@altsplices_seen, "retained intron");
#    }
#}

if(($add_string =~ /0+2+(1+)2+0+/) || ($add_string =~ /0+2+(1+)2+1*$/) || ($add_string =~ /^1*2+(1+)2+0+/) || ($add_string =~ /^1*2+(1+)2+1*$/)) {
    if(length($1) > 10) {
	$space_flag = 1;
	$report_string = $gtf_file1."\t"."retained intron\t\t\t";
	push(@altsplices_seen, "retained intron");
    }
}

if(@altsplices_seen != 0) {
    
    $report_string = join(":", @altsplices_seen);
    
#foreach my $thing (@altsplices_seen) {
#    print $thing, "\n";
#}
    
    print $gtf_file1, "\t", $report_string;
    
}
