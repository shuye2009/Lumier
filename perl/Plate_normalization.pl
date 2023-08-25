#!/usr/bin/perl -w

#Author: Shuye Pu
#Created on: Oct 18, 2016
#last modified on:

#working directory C:\RSYNC\LUMIER_C2H2\NOV232016

use Statistics::Basic qw(:all);
use Regexp::Common;

my @preys = ("KDM5C","ZKSCAN2", "ZNF76", "ZNF143", "ZNF238", "ZNF250", "ZNF554", "ZNF669", 
			"ZNF786", "ZNF784", "ZNF768", "ZNF16", "C13ORF8", "JARID1B", "ZNF512", "ZNF777", "ZNF704", 
			"ZNF45", "ZBTB42", "ZBTB26", "ZNF71", "ZNF192", "ZNF273", "ZNF317", "ZNF512_2",
			"ZNF768_2", "EGR2", "YY1", "ZNF112", "ZNF99", "ZNF280B", "ZNF281", "ZNF324B", "ZNF502",
			"ZNF574", "ZNF816A", "ZNF837", "JMJD5", "ZNF157", "ZNF280A"); # "PRDM2", 
@preys = sort(@preys);
my @treatments = ("NT", "RNASE"); #,"Rnase", "Notreatment"


foreach my $treat (@treatments){
	my %prey_bait_avg = ();
	my %prey_bait_stdev = ();
	my %baits = ();
	my $r_file = $treat."_replicates_and_average_raw_lum.tab";
	open(R, ">$r_file");
	print R "prey\tbait\tmean_lum\tstdev_lum\treplicates\n"; 
	foreach my $prey (sort @preys){
		my $input = $prey."_".$treat.".tab";
		
		my %p384_wellid_gene = (); #p384_wellid as key, gene as value
		my %gene_normalized_lum = (); #gene as key, array of normalized lums as value
		my %gene_AVGlum = (); #gene as key, average of lums as value
		my %p384_wellid_p96_lum = ();	#p384_wellid as key, lum as value, one key per well 
		my %p384_wellid_p96_lum_normalized = ();	#p384_wellid as key, lum as value, one key per well 
		my %p384_p96_lum = (); #p384_p96 as key, array of lums as value, one key per 96-well plate
		my %p384_p96_median = (); #p384_p96 as key, median of lums as value, , one key per 96-well plate
		
		#print "processing $input\n";	
		open(IN, "<$input") or die "can not open $input";
		<IN>;
		<IN>;
		while(<IN>){
			chomp($_);
		    my $temp = $_;
		    $temp =~ s/\n//;
		    $temp =~ s/\r//;
		
			
		    my @tary = split(/\t/,$temp);
		    
		    my $size = scalar(@tary);
		    
		    #if($input eq "JARID1B_RNASE.tab"){print "$input\t$size\t$temp\n";}
		    unless(defined $tary[0]){next;}
			if($tary[0] eq "" or $tary[0] eq "Plate Number"){next;}
			
			my $plate384 = $tary[0];
			my $wellID = $tary[1];
			my $well384 = $tary[2];
			my $lum = $tary[3];
			my $plate96 = $tary[4];
			my $well96 = $tary[5];
			my $sort_rep = $tary[6];
			my $gene = $tary[7];
			my $construct = $tary[8];
			
			if(defined $construct){
				$gene = $gene."_".$construct;
			}
			
			
			unless(defined $gene){$gene = "Blank";}
			
			if($lum eq "OVRFLW"){
				$lum = 5000000;
			}elsif($lum eq "X"){
				$lum = 1;
			}
			$baits{$gene} = 1;
			
			$p384_wellid_gene{"$plate384\t$wellID"} = $gene;
			$p384_wellid_p96_lum{"$plate384\t$wellID\t$plate96"} = $lum;
			#if($lum =~ m{^$RE{num}{int}$}){
			push(@{$p384_p96_lum{"$plate384\t$plate96"}}, $lum);
			#}
			
			
			
			
			
		}
		
		close(IN);
		
		while(my ($key, $ref) = each %p384_p96_lum){
			my @ar = @$ref;
			#my $half_size = int(scalar(@ar)/2);
			my @arr = sort{$a<=>$b}@ar;
			#my $median = $arr[$half_size];
			my $median = median($ref);
			$median =~ s/,//g;
			$p384_p96_median{$key} = $median;
			#print "$input\t$key\t$median\n";
		}
		
		my $ind = 0;
		while(my ($p384_wellid_p96, $lum) = each %p384_wellid_p96_lum){
			
			my @ids = split(/\t/, $p384_wellid_p96);
			my $p384 = $ids[0];
			my $wellid = $ids[1];
			my $p96 = $ids[2];
			my $median = $p384_p96_median{"$p384\t$p96"};
			unless(defined $median){next;}
			my $normalized_lum = 0;
			if($median ne 0){
				$normalized_lum = $lum/$median;
			}
			$normalized_lum = $lum; #output raw reading
			$p384_wellid_p96_lum_normalized{$p384_wellid_p96} = $normalized_lum;
			my $gene = $p384_wellid_gene{"$p384\t$wellid"};
				
			push(@{$gene_normalized_lum{$gene}}, $normalized_lum);
			
			$ind++;
			#print "$ind\t$treat\t$prey\t$p384\t$wellid\t$p96\t$gene\t$lum\t$median\t$normalized_lum\n";
			
		}
		
		foreach my $gene (sort(keys %gene_normalized_lum)){
			my $aryref = $gene_normalized_lum{$gene};
			
			my $mean_lum = mean($aryref);
			print "before $mean_lum\n";
			$mean_lum =~ s/,//g;
			print "after $mean_lum\n";
			my $stdev_lum = stddev($aryref);
			$prey_bait_avg{$prey}->{$gene} = $mean_lum;
			$prey_bait_stdev{$prey}->{$gene} = $stdev_lum;
			my @l = @$aryref;
			my $sl = join(", ", @l);
			print R "$prey\t$gene\t$mean_lum\t$stdev_lum\t$sl\n";
		}
	}
	
	close(R);
	
	my $out = $treat."_bait_prey_raw_matrix.tab";
	my $out1 = $treat."_bait_prey_raw_matrix_cluster.tab";
	my @all_baits = sort(keys %baits);
	
	open(OUT, ">$out");
	open(OUT1, ">$out1");
	my $p = join("\tSD\t", @preys);
	my $p1 = join("\t", @preys);
	print OUT "\t$p\n";
	print OUT1 "\t$p1\n";
	foreach my $bait(@all_baits){
		print OUT "$bait";
		print OUT1 "$bait";
		
		foreach my $prey (@preys){
			my $v = "NaN";
			my $sd = "NaN";
			if(exists $prey_bait_avg{$prey}->{$bait}){
				$v = $prey_bait_avg{$prey}->{$bait};	
			}
			if(exists $prey_bait_stdev{$prey}->{$bait}){
				$sd = $prey_bait_stdev{$prey}->{$bait};	
			}
			
			print OUT "\t$v\t$sd";
			print OUT1 "\t$v";
		}
		print OUT "\n";
		print OUT1 "\n";
	}
	close(OUT);
	close(OUT1);
}



print "finish";