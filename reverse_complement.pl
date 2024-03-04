#! /usr/bin/perl -w
use strict;

# hash used to translate DNA symbols into complements
my %tr_hash = qw (
		  A T
		  a t
		  T A
		  t a
		  G C
		  g c
		  C G
		  c g
		  N N
		  n n
		  M K
		  m k
		  R Y
		  r r
		  W W
		  w w
		  S S
		  s s
		  K R
		  k r
		  Y M
		  y m
		  );


my($filename) = @ARGV;

my $usage = "USAGE:  ./reverse_complement <fasta_file>
reverse complements sequences in fasta or multi 
fasta file containing DNA symbols, 2 place 
ambiguity codes and Ns (everything but 3 place 
ambiguity codes)\n";

die $usage unless (@ARGV == 1);

my %hash;
read_fasta($filename, \%hash);

# note:  this will not maintain the ordering of the original multi fasta file
foreach my $header (keys %hash) {

    my $sequence = $hash{$header};

    my $reverse_sequence = reverse($sequence);

    my $reverse_complement = "";

    my $length = length($sequence);

    # yeah it's inefficent but it gets the job done
    for(my $i=0;$i<$length;$i++) {
	my $character = substr($reverse_sequence, $i, 1);

	if(!defined($tr_hash{$character})) {
	    die("unexpected DNA symbol encountered $character... exiting");
	}

	$reverse_complement .= $tr_hash{$character};

    }

    print ">${header}_reverse_complement\n";
    print $reverse_complement, "\n";

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
