use ExtUtils::MakeMaker;
# See lib/ExtUtils/MakeMaker.pm for details of how to influence
# the contents of the Makefile that is written.
WriteMakefile(
    NAME              => 'trace_recalling',
    VERSION      => '0.5',
    PREREQ_PM         => {},
    EXE_FILES      => [ qw(est_genome_ambig compare_altsplices_single.pl est2geno_to_gtf.pl trace_recalling.pl reverse_complement.pl) ],
);  
