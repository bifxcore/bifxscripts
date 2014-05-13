#!/depot/perl-5.12.1/bin/perl
#
# PURPOSE: To create a SL maps from converted GFF table and converted .alignTable 
# 
# USAGE: $0 <Gff Table(from sl2_map.pl)> <Align Table(from sl3_map.pl)> <outputlevel[gene|slsite]>
# 
# Output: will be in <STDOUT>. Ignore Warnings
#
###############################################################################
use strict;
use strict;
use warnings ;
use Pod::Usage;
use Data::Dumper;

#use PerlIO::gzip;



#-------------------------------------------------------------------------------
# Main.
#-------------------------------------------------------------------------------
#my $gffExtra='apidb\|Lmjchr';
#my $alignExtra ='psu\|Lmjchr';
#my $feture2Capture ='gene';
my $sumStat;


my $gffFile = shift;
my $alignFile = shift;
my $outputType =shift;
if(!$outputType){$outputType = "gene";}

my (%gff, %byId, %chromosome) = (); 
open(READ, $gffFile);# die "Could not read $gffFile File";
while(my $gffLine = <READ>){
	chomp($gffLine);
	my ($gffId, $gffChr, $gffStrand, $gffUtrSt, $gffCdsStart, $gffCdsEnd) = split("\t", $gffLine);
	$gff{$gffChr}{$gffId}=[$gffId, $gffStrand, $gffUtrSt, $gffCdsStart, $gffCdsEnd];
	$byId{$gffId}=[$gffChr, $gffStrand, $gffUtrSt, $gffCdsStart, $gffCdsEnd];
	unless(exists $chromosome{$gffChr}){
		$chromosome{$gffChr}="1";
	}
}
close(READ);

my (%result, %antiresult, %alignHash, %ranks) = ();
open(PADI, $alignFile);
while(my $alignLine= <PADI>){
	chomp($alignLine);
	#print "$alignLine\n";
    my($alignChr, $alignPos, $posAlignCount, $negAlignCount) = split ("\t", $alignLine);
    #store in a hash for later use (to count whole chro hits)
    $alignHash{$alignChr}{$alignPos}=$posAlignCount;
    $alignHash{$alignChr}{$alignPos}=$negAlignCount;
    
    #search
    
    for my $geneId (keys %{$gff{$alignChr}}){
    	my $utrStart=$gff{$alignChr}{$geneId}[2];
    	my $cdsStart=$gff{$alignChr}{$geneId}[3];
    	my $cdsEnd  =$gff{$alignChr}{$geneId}[4];
    	if($alignPos > $utrStart && $alignPos <= $cdsEnd){
    		$result{$geneId}{$alignPos}=$posAlignCount if $posAlignCount > 0;
    		$antiresult{$geneId}{$alignPos}=$negAlignCount if $negAlignCount > 0;
    	}
    	if($alignPos >= $cdsStart && $alignPos < $utrStart){
    		$result{$geneId}{$alignPos}=$negAlignCount if $negAlignCount > 0;
    		$antiresult{$geneId}{$alignPos}=$posAlignCount if $posAlignCount > 0;
    	}
    }
}

#Find Total reads per gene/chr
my $hrefTotalPerGene= &countHitPerGene(\%result);
my $hrefTotalPerChr = &countTotalPerChr(\%alignHash);
&find_ranks();

my $hrefResultsWithName=&nameSLsites02(\%result);
if($outputType =~ /slsite/){&printResultsSLsites($hrefResultsWithName, $hrefTotalPerChr, $hrefTotalPerGene);}
if($outputType =~ /gene/){&printResultsGene($hrefResultsWithName, $hrefTotalPerChr, $hrefTotalPerGene);}

#-------------------------------------------------------------------------------
# Functions
#-------------------------------------------------------------------------------

sub countHitPerGene($){
    my $href = shift;
    my %totalHitsPerGene = ();
    for my $geneId (sort keys %$href){
        for my $pos (sort {$a<=>$b} keys %{$href->{$geneId}}){
            $totalHitsPerGene{$geneId} += $$href{$geneId}{$pos};
        }
    }
    return \%totalHitsPerGene;	
}

sub countTotalPerChr($){
    my $href = shift;
    my %totalHitsPerChr = ();
    for my $chr (sort {$a<=>$b} keys %$href){
	for my $pos (keys %{$href -> {$chr}}){
	    $totalHitsPerChr{$chr} += $$href{$chr}{$pos};
	}
    }
    return \%totalHitsPerChr;
}

sub nameSLsites02(){
    my %newResults =();
    for my $geneId (keys %result){
	my($intervel); my (@upInts, @downInts);
	
	# for genes in + strand
	if($byId{$geneId}[1] eq "+"){
	    #gene is on positive strand of ref
	    for my $pos (sort {$a<=>$b} keys %{$result{$geneId}}){
			my $realCDSstart = $byId{$geneId}[3];
			$intervel = $pos - $realCDSstart;	
			if($intervel < 0){
			    #alignPos is upstream to ATG(cdsStart)
		    	unshift (@upInts, $intervel);
			    #newResults{geneids}{position}=[geneStrand, alignPos alignCount, chr, Name(to be assigned) ]
			    $newResults{$geneId}{$intervel}=[$byId{$geneId}[1], $pos, $result{$geneId}{$pos}, $byId{$geneId}[0], $ranks{$geneId}{$pos}];
			}else{
		    	#alignPos is dowstream to ATG(cdsStart)
			    push(@downInts, $intervel);
			    $newResults{$geneId}{$intervel}=[$byId{$geneId}[1], $pos, $result{$geneId}{$pos}, $byId{$geneId}[0], $ranks{$geneId}{$pos}];
			}
	    }
	    
	    # done accumulating data. Now name it
	    my $maxUpPoints=@upInts; my $pointCount =0;
	    foreach my $upInt (@upInts){
			++$pointCount;
			my $pointName = "SLUp_$pointCount";
			$newResults{$geneId}{$upInt}[5]=$pointName;
			#print "$geneId\t$newResults{$geneId}{$upInt}[1]\t$pointName\n";
	    }
	    my $maxDownPoints=@downInts; $pointCount =0;
        foreach my $downInt (@downInts){
            ++$pointCount;
            my $pointName = "SLDo_$pointCount";
            $newResults{$geneId}{$downInt}[5]=$pointName;
            #print "$geneId\t$newResults{$geneId}{$downInt}[1]\t$pointName\n";
        }
	}else{
	    # for genes in - strand
	    for my $pos (sort {$a<=>$b} keys %{$result{$geneId}}){
			my $realCDSstart = $byId{$geneId}[4];
			$intervel = $realCDSstart - $pos;
			if($intervel > 0){
			    #alignPos is downstream to ATG(cdsStart)
			    unshift(@downInts, $intervel);
			    $newResults{$geneId}{$intervel}=[$byId{$geneId}[1], $pos, $result{$geneId}{$pos}, $byId{$geneId}[0], $ranks{$geneId}{$pos}];
			}else{
			 	#alignPos is upstream to ATG(cdsStart)
			    push(@upInts, $intervel);
			    $newResults{$geneId}{$intervel}=[$byId{$geneId}[1], $pos, $result{$geneId}{$pos}, $byId{$geneId}[0], $ranks{$geneId}{$pos}];
			}  	
	    }
	    # done accumulating data. Now name it
        my $maxUpPoints=@upInts; my $pointCount =0;
        foreach my $upInt (@upInts){
            ++$pointCount;
            my $pointName = "SLUp_$pointCount";
            $newResults{$geneId}{$upInt}[5]=$pointName;
            #print "$geneId\t$newResults{$geneId}{$upInt}[1]\t$pointName\n";
        }
        my $maxDownPoints=@downInts; $pointCount =0;
        foreach my $downInt (@downInts){
            ++$pointCount;
            my $pointName = "SLDo_$pointCount";
            $newResults{$geneId}{$downInt}[5]=$pointName;
            #print "$geneId\t$newResults{$geneId}{$downInt}[1]\t$pointName\n";            
        }
	}
	
	}

	
	
    
    return \%newResults;
}


sub find_ranks(){
	for my $geneId ( sort keys %result){	
		my $rank = 1; my $prevCount = "0";
	    for my $pos (sort {${$result{$geneId}}{$b}<=>${$result{$geneId}}{$a}} keys %{$result{$geneId}}){
	    	my $count = $result{$geneId}{$pos};
	    	my $iRank = $rank;
	    	if ($count == $prevCount){
	    		--$rank;
	    		$iRank = $rank.".1";	
	    	}
	    	$prevCount = $count;
	    	$ranks{$geneId}{$pos}=$iRank;
			++$rank;
	   }
	}
	return;
}


sub printResultsSLsites(){
    #newResult{geneids}{intervels}=[geneStrand, alignPos, alignCount, chr, Name(to be assigned) ]
    print "#SL ID\tChr\tAlign Position\tAlign Count\tSENSE\tGene ID\tStrand\tCDSstart\tCDSend\tPutative 5' end\tSL Rank\tIntervel\tPeak percent\n";
    my $href1 = shift;
    my $hrefChrCount =shift;
    my $hrefGeneCount = shift;
    for my $geneId (sort keys %$href1){
        #set correct CDS start:
    	my $geneStrand = $byId{$geneId}[1];
	    my $noof_SLsitesPerGene= keys %{$href1 ->{$geneId}};
    	my $realCDSstart; my $realCDSend; my $realPrevGeneEnd;
	    if($geneStrand eq "+"){
	    	$realCDSstart = $byId{$geneId}[3];
	    	$realCDSend = $byId{$geneId}[4];
	    	$realPrevGeneEnd = $byId{$geneId}[2];
	    }
	    else{
	       $realCDSstart = $byId{$geneId}[4];
	       $realCDSend = $byId{$geneId}[3];
	       $realPrevGeneEnd = $byId{$geneId}[2];
	    }

	    for my $int (sort {$a<=>$b} keys %{$href1 ->{$geneId}}){
	       #my $geneStrand = $$href1{$geneId}{$int}[0];
	       my $alignPos = $$href1{$geneId}{$int}[1];
    	   my $alignCount = $$href1{$geneId}{$int}[2];
	       my $peakRatio = (($alignCount*100)/$$hrefGeneCount{$geneId});
	       my $alignChr = $$href1{$geneId}{$int}[3];
    	   my $slRank =  $$href1{$geneId}{$int}[4];
	       my $slName = $$href1{$geneId}{$int}[5];
	       my $slId = $geneId."_".$slName;
	       my $strand = 1; $strand = "-1" if $geneStrand eq '-';
           print "$slId\t$alignChr\t$alignPos\t$alignCount\tSENSE\t$geneId\t$strand\t$realCDSstart\t$realCDSend\t$realPrevGeneEnd\t$slRank\t$int\t$peakRatio\n";
	    }
	    # Now try printing antisense reads for that gene
	    foreach my $antipos (keys %{$antiresult{$geneId}}){
	    	my $strand = 1; $strand = "-1" if $geneStrand eq '-';
	    	my $int = $realCDSstart - $antipos;
	    	my $alignChr = $byId{$geneId}[0];
	    	my $slid = $geneId."_antisl";
	    	print "$slid\t$alignChr\t$antipos\t$antiresult{$geneId}{$antipos}\tANTI\t$geneId\t$strand\t$realCDSstart\t$realCDSend\t$realPrevGeneEnd\t0\t$int\t0\n";
	    }
	    
    }#for my $geneId
}

sub printResultsGene(){
	print "#GeneId\tStrand\tChr\tCDS start\tCDS end\tPutative 5' End\t#SLsitesPerGene\tAntiSiteCount\tSenseReadCount\tAntiSiteReadcount\t#upSites\t#downSites\n";
    #newResult{geneids}{intervels}=[geneStrand, alignPos alignCount, chr, Name(to be assigned) ]
    my $href1 = shift;
    my $hrefChrCount =shift;
    my $hrefGeneCount = shift;
    my (@allMajorPeakCounts,@allMajorPeakRatios,@allMajorPeakPos_fromATG,@allSecondMajorPeakCounts,@allSecondMajorPeakRatios,@allSecondMajorPeakPos_fromATG, @allPeakRatioDiffs);
    my %bySlId=();
    
    for my $geneId (sort keys %$href1){
	#set correct CDS start:
	my $geneStrand = $byId{$geneId}[1];
	my $noof_SLsitesPerGene= keys %{$href1 ->{$geneId}};
	my ($realCDSstart, $realCDSend, $realPrevGeneEnd);
	my $noof_upSites = "0"; my $noof_doSites = "0";
	if($geneStrand eq "+"){
            $realCDSstart = $byId{$geneId}[3];
            $realCDSend = $byId{$geneId}[4];
            $realPrevGeneEnd = $byId{$geneId}[2];
    }else{
           $realCDSstart = $byId{$geneId}[4];
           $realCDSend = $byId{$geneId}[3];
           $realPrevGeneEnd = $byId{$geneId}[2];
    }
	my $alignChr = $byId{$geneId}[0];

	my %sl_site_info=();
	for my $int (sort {$a<=>$b} keys %{$href1 ->{$geneId}}){
	    my $alignPos = $$href1{$geneId}{$int}[1];
	    my $alignCount = $$href1{$geneId}{$int}[2];
	    my $peakRatio = (($alignCount*100)/$$hrefGeneCount{$geneId});
	    my $slName = $$href1{$geneId}{$int}[5];
	    my $slId = $geneId."_".$slName;
	    
	    #find total up/down sl site info
	    if($int < 0){++$noof_upSites;}else{++$noof_doSites;}
	    #find first up/down sl site info
	    my $upslid_exp="SLUp_1";my $doslid_exp="SLDo_1";
	    if($upslid_exp eq $slName){
	    	$bySlId{"firstUpId"}=$slId;
	    	$bySlId{"firstUpPos"}=$int;
	    	$bySlId{"firstUpCount"}=$alignCount;
	    	$bySlId{"firstUpRatio"}=$peakRatio;
	    	
	    }
	    if($doslid_exp eq $slName){
	    	$bySlId{"firstDoId"}=$slId;
	    	$bySlId{"firstDoPos"}=$int;
	    	$bySlId{"firstDoCount"}=$alignCount;
	    	$bySlId{"firstDoRatio"}=$peakRatio;
	    }
	    
	    $sl_site_info{$peakRatio}=[$slName, $slId, $alignPos, $alignCount];
	}
	my @peakRatios = sort{$b<=>$a} keys %sl_site_info;
	my $majorPeakId = $sl_site_info{$peakRatios[0]}[1];
	my $majorPeakPos_fromATG;
	 if($geneStrand eq "+"){
	 	$majorPeakPos_fromATG = $sl_site_info{$peakRatios[0]}[2] - $realCDSstart ;
	 }else{
	    $majorPeakPos_fromATG = $realCDSstart - $sl_site_info{$peakRatios[0]}[2];
	 }
	my $majorPeakCounts = $sl_site_info{$peakRatios[0]}[3];      
        my $majorPeakRatio = (($majorPeakCounts*100) / $$hrefGeneCount{$geneId});
        
        my($secondMajorPeakId,$secondMajorPeakPos_fromATG,$secondMajorPeakCounts,$secondMajorPeakRatio,$peakRatioDiff);
        if($noof_SLsitesPerGene > "1"){
            $secondMajorPeakId = $sl_site_info{$peakRatios[1]}[1];
	   		$secondMajorPeakPos_fromATG = $realCDSstart - $sl_site_info{$peakRatios[1]}[2];
	    	$secondMajorPeakCounts = $sl_site_info{$peakRatios[1]}[3];      
            $secondMajorPeakRatio = (($secondMajorPeakCounts*100) / $$hrefGeneCount{$geneId});
            $peakRatioDiff = $majorPeakRatio - $secondMajorPeakRatio;
        }
        #okay now find out antisense reads if any (from %antiresult)
        my $antiSiteCount = "0";my $antiSite_readcount = "0";
        foreach my $antipos (keys %{$antiresult{$geneId}}){
        	++$antiSiteCount;
        	$antiSite_readcount += $antiresult{$geneId}{$antipos};
        }
        #print 
        my $strand = 1;
        $strand = "-1" if $geneStrand eq '-';
		print "$geneId\t$strand\t$alignChr\t$realCDSstart\t$realCDSend\t$realPrevGeneEnd\t$noof_SLsitesPerGene\t$antiSiteCount\t$$hrefGeneCount{$geneId}\t$antiSite_readcount\t$noof_upSites\t$noof_doSites\n";
		
		
       }#for my $geneId
  
  #do summary statistics (obsolete)
  if($sumStat){
    printf("MajPkCount Ave =%d\nMajPkCount median = %d\n",  ave(\@allMajorPeakCounts), median(\@allMajorPeakCounts));
    printf("MajPkRatio Ave =%d\nMajPkRatio median = %d\n", ave(\@allMajorPeakRatios), median(\@allMajorPeakRatios));
    printf("MajPkPos Ave =%d\nMajPkPos median = %d\n", ave(\@allMajorPeakPos_fromATG), median(\@allMajorPeakPos_fromATG));
  
    printf("SecMajPkCount Ave =%d\nSecMajPkCount median = %d\n",  ave(\@allSecondMajorPeakCounts), median(\@allSecondMajorPeakCounts));
    printf("SecMajPkRatio Ave =%d\nSecMajPkRatio median = %d\n", ave(\@allSecondMajorPeakRatios), median(\@allSecondMajorPeakRatios));
    printf("SecMajPkPos Ave =%d\nSecMajPkPos median = %d\n", ave(\@allSecondMajorPeakPos_fromATG), median(\@allSecondMajorPeakPos_fromATG));
  
    printf("PkRatioDiff Ave =%d\nPkRatioDiff median=%d", ave(\@allPeakRatioDiffs), meidan(\@allPeakRatioDiffs));
    
  }
}

sub median(){
	my $aRef = shift;
	my $median;
	my @sortedArray = sort{$a<=>$b} @$aRef;
	my $aSize = $#sortedArray+1;
    open(my $fh, ">>array");
	foreach (@sortedArray){print $fh "$_\n";}
	if ($aSize % 2){
		$median = $sortedArray[int($aSize/2)];
	}else{
		$median = ($sortedArray[$aSize/2] + $sortedArray[$aSize/2-1]) / 2;
	}
	return $median;
}

sub ave(){
	my $aRef =shift;
	my $sum = "";
	my $aSize = @$aRef;
	foreach (@$aRef){
		$sum += $_;
	}
	my $ave = $sum/$aSize;
	return $ave;
}
