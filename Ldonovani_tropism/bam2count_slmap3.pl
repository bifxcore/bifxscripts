#!/depot/perl-5.12.1/bin/perl
#!/depot/ergatis-perl-5.8.9/bin/perl
use strict;
use warnings;
use Data::Dumper;
use Bio::DB::Sam;
#use Array::IntSpan;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
#----------------------------------------------------------------------
# handle command line
#----------------------------------------------------------------------
my %options =();
my $opt_obj = GetOptions (\%options,
              'bam|b=s',
              'help|h',
              'log|l=s') || pod2usage();

my $usage ='$0 -bam=<bam file>';
print "$usage\n" if($options{help});

my $debug = 0;$debug = 2 if($options{debug});

#----------------------------------------------------------------------
# Main - Read input
#----------------------------------------------------------------------


# Now count pileups in  bam
my $sam = Bio::DB::Sam->new(-bam => $options{bam});
my @contigs = $sam->seq_ids;
foreach my $contig (@contigs){
    my %counts= ();
	my @alignments = $sam->get_features_by_location(
                    -seq_id  => $contig);
    for my $a (@alignments){
    	my $seqid = $a->seq_id;
        my $start = $a->start;
        my $end   = $a->end;
        my $strand= $a->strand;
        my $posi;
        $posi = $start if $strand > 0;
        $posi = $end if $strand < 0;
       if (exists $counts{$posi}{$strand}){
       	++$counts{$posi}{$strand};
       }else{
       	$counts{$posi}{$strand}="1";
       }
    }
    #print Dumper(\%counts);
    foreach my $pos (sort {$a <=> $b} keys %counts){
        my $pos_count = 0;
        my $neg_count = 0;
        $pos_count = $counts{$pos}{'1'} if (exists $counts{$pos}{'1'});
        $neg_count = $counts{$pos}{'-1'} if (exists $counts{$pos}{'-1'});
        #$contig =~ s/\D//g;
    	print "$contig\t$pos\t$pos_count\t$neg_count\n";
    }
}
