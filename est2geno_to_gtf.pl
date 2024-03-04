#! /usr/bin/perl -w

use strict;
use EGlite;
use FAlite;
use POSIX;

my($filename) = @ARGV;

#my$dir = $ARGV[0];
#mask_genoseq($dir);



#sub mask_genoseq{



#  my($tmpdir) = @_;
my@mask_txt = ();

# read in all *.out2 files (those contain my est2geno alignments, I guess you should change this)
#  opendir (DIR, "$tmpdir");
#  while (defined(my$file = readdir(DIR))){

#print "new file here\n";

#    unless($file =~ /.*\.out2$/){next}

open(INF, $filename);

#open (INF, "$tmpdir\/$file")||die "cannot open file";

if(my $report = new EGlite(\*INF)){
    my$genome = $report->genome;
    my($est) = $report->est;
    my$orientation= '-';
    my$orient = $report->genomeDirection;
    if($orient eq "forward"){
        $orientation = '+';
    }
    my%maskexons = ();
    my@intron_sizes=();
    my@exon_sizes=();
    
    while(my $element = $report->nextElement) {
# skip introns
        if($element->type eq "INTRON"){
	    next;
        }
	
# the other elements are exons
        my$start = $element->genomeBegin;
        my$end = $element->genomeEnd;
	
	print "sometext\tPMASK\tCDS\t$start\t$end\t\.\t$orientation\t\.\tgene_id \"${filename}\"\; transcript_id \"${filename}\"\;\n";
	
        push (@mask_txt, 
	      "sometext\tPMASK\tCDS\t$start\t$end\t\.\t$orientation\t\.\tgene_id \"${filename}\"\; transcript_id \"${filename}\"\;\n");
    }
    
# if there is no readable file in the *out2: issue warning
}else{
#      print STDERR "$tmpdir\/$file is empty or has wrong format\n"; 
}
close INF;
#  }
#  closedir DIR;

#  print @mask_txt;
#}


