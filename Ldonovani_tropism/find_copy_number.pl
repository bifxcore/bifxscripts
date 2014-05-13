#!/depot/perl-5.12.1/bin/perl
#
#  Calculate coverage (median/average) from BAM file
#  Usage $0 -bam=<bam file>
#
#
##############

use strict;
use warnings;
use Bio::DB::Sam;
use Math::NumberCruncher;
use Data::Dumper;
use Math::Round;

use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);


#----------------------------------------------------------------------
# handle command line
#----------------------------------------------------------------------
my %options =();
my $opt_obj = GetOptions (\%options,
              'bam|b=s',
              'bin_size|s=s',
              'cov_type=s',
              'contig=s',
              'ploidy',
              'intervaldata=s',
              'intervaldata_withid=s',
              'debug|d',
              'help|h');

my $usage = qq{
$0 
	-b, --bam
		BAMFILE
		
	-s, --bin_size
		# of bases
		
	-cov_type
		median|mean|mode|null
		use null option to find out bases/regions with NO coverage at all.
		When null option is used, bin_size defaults to "1"
		And only bases with NO coverage will be reported.
		
	-contig
		Name of the contig. Default all contigs
		
		
	-ploidy
		Set this if polidy needs to be calculated on the fly
		Formulae ((cov_of_region/cov_for_entire_contig)*2)
		Should be set only with -bin_size
		Will not work for whole chromosomes.
	-intervaldata
	   A file with three tab delimited columns
	   <col1> : Contig name
	   <col2> : Start position
	   <col3> : End position
		
	-h, --help
	-d, --debug
};

my $example = qq{
	$0 -bam ES017r.sorted.bam -cov_type median > LindoOnly_vs_LinJ.cov.median
	   - gives one median coverage value per contig
	   - considers coverage at every base to calculate median
	
	$0 -bam ES017r.sorted.bam -cov_type median -s 1
	   - prints coverage per base for each position 
	   - outputs covearge of all contigs in bam unless specified by -contig tag
	   - good for drawing covearge plots. eg in user plots in artemis
	   
	grep ">" /core/bifx/NGS/bowtie2_index/LmjF-TritrypV3_3.fasta |perl -pi -e "s/>//" |xargs -i echo " $0 "-bam ES017r.sorted.bam -cov_type median -s 1000 -contig \{\} \> \{\}.covearge.plot
	   - same as above, just automates creating one file for each contig.
};

die("$usage\n") if $options{help};
die(print "$usage\n") if(!$options{bam});
#die("To use -polidy, set a -bin_size as well\n $usage") if ($options{ploidy} && !$options{bin_size});
$options{bin_size}="1"  if $options{cov_type} eq "null";
print "WARNING: Bin size set to 1 'cause cov_type requested is \"null\"\n";
#----------------------------------------------------------------------
# main
#----------------------------------------------------------------------

my $sam = Bio::DB::Sam->new(-bam => $options{bam});

# get all contig ids from BAM
my @contigs = $sam->seq_ids;

# if polidy is asked, get over all chromosome coverage first
my %contig_level_cov = ();my $ploidy ;
if($options{ploidy}){
	foreach my $contig (@contigs){
		next if ($options{contig} && $options{contig} ne $contig);
		my $cov_value ;
		my $contiglen = $sam->length($contig);
		my $start =1;
		$cov_value = &find_coverage($sam, $contig, $start, $contiglen, $options{cov_type});
		$contig_level_cov{$contig}=$cov_value;	
	}	
}
print Dumper(\%contig_level_cov) if $options{debug};


# if contigname, intervel are supplied from a file
if($options{intervaldata}){
    # read interval file:
    open(READ, $options{intervaldata});
    my %intdata = ();
    my $hashkey = "1";
    while(my $line = <READ>){
    	chomp($line);
    	#print $line;
    	my ($contig, $start, $end)= split("\t", $line);
        $intdata{$hashkey}{'contig'}=$contig;
        $intdata{$hashkey}{'start'}=$start;
        $intdata{$hashkey}{'end'}=$end;
        ++$hashkey;
    }
	my %results=();
	foreach my $key (keys %intdata){
	    my $contig = $intdata{$key}{contig};
	    # skip unwanted contigs if one asked specifically
	    next if ($options{contig} && $options{contig} ne $contig);
	
	    my $cov_value;
	    # Fix segment length from set bin_size
	    my $start = $intdata{$key}{start}; 
	    my $end = $intdata{$key}{end};
	    my $contiglen = $sam->length($contig);
	    #print "Global:$contig\t$start\t$end\n";

	    if(!$options{bin_size}){
	        $cov_value = &find_coverage($sam, $contig, $start, $end, $options{cov_type});
	        print "$contig\t$start\t$end\t$cov_value\n";
	    }else{
	    	my $oristart = $start;
	    	my $oriend = $end;
	        for(my $x = $oristart ; $x < $oriend; $x = $x + $options{bin_size}){
	           $start = $x+1;
	           next if $start >= $contiglen;
	           $end = $start + $options{bin_size} - 1;
	           $end = $oriend if ($end > $oriend);
	           $end = $contiglen if ($end > $contiglen);
	           #print "Actual2:$contig\t$start\t$end\n";
	            
	            $cov_value = &find_coverage($sam, $contig, $start, $end, $options{cov_type});
	            $ploidy = ($cov_value/$contig_level_cov{$contig})*2 if $options{ploidy};
	            print "$contig\t$start\t$end\t$cov_value\n" if (!$options{ploidy} && $options{cov_type} ne "null");
	            print "$contig\t$start\t$end\t$cov_value\n" if ($options{cov_type} eq "null" && $cov_value == "0");
	            print "$contig\t$start\t$end\t$cov_value\t$ploidy\n" if $options{ploidy};
	            
	            $start = $end+1;
	            $end += $options{bin_size};
	        }
	    }
	}
}elsif($options{intervaldata_withid}){
 # if contigname, interveldata with id are supplied from a file
 
    # read interval file:
    open(READ, $options{intervaldata_withid});
    my %intdata = ();
    my $hashkey = "1";
    while(my $line = <READ>){
        chomp($line);
        #print $line;
        my ($id, $contig, $start, $end)= split("\t", $line);
        $intdata{$hashkey}{'id'}=$id;
        $intdata{$hashkey}{'contig'}=$contig;
        $intdata{$hashkey}{'start'}=$start;
        $intdata{$hashkey}{'end'}=$end;
        ++$hashkey;
    }
    my %results=();
    foreach my $key (keys %intdata){
        my $contig = $intdata{$key}{contig};
        # skip unwanted contigs if one asked specifically
        next if ($options{contig} && $options{contig} ne $contig);
    
        my $cov_value;
        # Fix segment length from set bin_size
        my $id = $intdata{$key}{id};
        my $start = $intdata{$key}{start}; 
        my $end = $intdata{$key}{end};
        my $contiglen = $sam->length($contig);
        $end = $contiglen if $end > $contiglen;
        print "Global:$id\t$contig\t$start\t$end\n";

        if(!$options{bin_size}){
            $cov_value = &find_coverage($sam, $contig, $start, $end, $options{cov_type});
            print "$id\t$contig\t$start\t$end\t$cov_value\n" if !$options{ploidy};
            $ploidy = ($cov_value/$contig_level_cov{$contig})*2 if $options{ploidy};
           
            print "$id\t$contig\t$start\t$end\t$cov_value\t$ploidy\n" if $options{ploidy};
            
        }else{
            my $oristart = $start;
            my $oriend = $end;
            for(my $x = $oristart ; $x < $oriend; $x = $x + $options{bin_size}){
               $start = $x+1;
               next if $start >= $contiglen;
               $end = $start + $options{bin_size} - 1;
               $end = $oriend if ($end > $oriend);
               $end = $contiglen if ($end > $contiglen);
               #print "Actual2:$contig\t$start\t$end\n";
                
                $cov_value = &find_coverage($sam, $contig, $start, $end, $options{cov_type});
                $ploidy = ($cov_value/$contig_level_cov{$contig})*2 if $options{ploidy};
                print "$contig\t$start\t$end\t$cov_value\n" if (!$options{ploidy} && $options{cov_type} ne "null");
                print "$contig\t$start\t$end\t$cov_value\n" if ($options{cov_type} eq "null" && $cov_value == "0");
                print "$contig\t$start\t$end\t$cov_value\t$ploidy\n" if $options{ploidy};
                
                $start = $end+1;
                $end += $options{bin_size};
            }
        }
    }
}else{
	# Default original code. This part works when no input file is given, i.e entire bam files is processed
	# loop to split each contig
	my %results=();
	foreach my $contig (@contigs){
	    # skip unwanted contigs if one asked specifically
	    next if ($options{contig} && $options{contig} ne $contig);
	
	    my $cov_value ;
	    #create segments of set bin size
	    # find length of entire segment
	    my $contiglen = $sam->length($contig);
	    
	    # Fix segment length from set bin_size
	    my $start; my $end;
	    $start = 1;
	    #print "$contig\t$contiglen\n";
	    if(!$options{bin_size}){
	        $end = $contiglen; #consider full len if bin_size is not set
	        $cov_value = &find_coverage($sam, $contig, $start, $end, $options{cov_type});
	        print "$contig\t$start\t$end\t$cov_value\n";
	    }else{
	        $end = $start + $options{bin_size} -1;
	        for(my $x = 0 ; $x < $contiglen; $x = $x + $options{bin_size}){
	            $end = $contiglen if ($end > $contiglen);
	            $cov_value = &find_coverage($sam, $contig, $start, $end, $options{cov_type});
	            $ploidy = ($cov_value/$contig_level_cov{$contig})*2 if $options{ploidy};
	            print "$contig\t$start\t$end\t$cov_value\n" if (!$options{ploidy} && $options{cov_type} ne "null");
	            print "$contig\t$start\t$end\t$cov_value\n" if ($options{cov_type} eq "null" && $cov_value == "0");
	            print "$contig\t$start\t$end\t$cov_value\t$ploidy\n" if $options{ploidy};
	            
	            $start = $end+1;
	            $end += $options{bin_size};
	        }
	    }
	}
		
}
#----------------------------------------------------------------------
# Functions
#----------------------------------------------------------------------
sub find_coverage(){
	my $sam = shift;
	my $contig = shift;
	my $start = shift;
	my $end = shift;
	my $cov_type = shift;

	my $segment = $sam->segment(-seq_id=>$contig, -start=>$start, -end=>$end);
	my ($coverage) = $segment->features('coverage');
    my @cov_data = $coverage->coverage;
    
    my $median = Math::NumberCruncher::Median(\@cov_data);
    my $mean = Math::NumberCruncher::Mean(\@cov_data);
    my $mode= Math::NumberCruncher::Mode(\@cov_data);
    my $rawtotal = "0";
    foreach my $val (@cov_data){
	$rawtotal += $val;
    }
    return ($median) if $cov_type eq "median";
    return ($mean) if $cov_type eq "mean";
    return ($mode) if $cov_type eq "mode";
    return ($rawtotal) if $cov_type eq "total";

    return $cov_data[0] if $cov_type eq "null";
}
           

