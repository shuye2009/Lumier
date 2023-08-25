#!/usr/bin/perl -w

#Author: Shuye Pu
#Created on: Oct 18, 2016
#last modified on: March 08, 2018

# working directory C:\RSYNC\LUMIER_C2H2\SecondProtocol
# This program implements the normalization methods detailed in:
# Methods Mol Biol. 2013;1064:17-27. doi: 10.1007/978-1-62703-601-6_2.
# Analysis of protein-protein interactions using LUMIER assays.
# Blasche S, Koegl M.

use Statistics::Descriptive qw(:all);
use Regexp::Common;



#my @preys = ("KDM5C","ZKSCAN2", "ZNF76", "ZNF143", "ZNF238", "ZNF250", "ZNF554", "ZNF669", 
#			"ZNF786", "ZNF784", "ZNF768", "ZNF16", "C13ORF8", "JARID1B", "ZNF512", "ZNF777", "ZNF704", 
#			"ZNF45", "ZBTB42", "ZBTB26", "ZNF71", "ZNF192", "ZNF273", "ZNF317", "ZNF512_2",
#			"ZNF768_2", "EGR2", "YY1", "ZNF112", "ZNF99", "ZNF280B", "ZNF281", "ZNF324B", "ZNF502",
#			"ZNF574", "ZNF816A", "ZNF837", "JMJD5", "ZNF157", "ZNF280A"); # "PRDM2", 

#			
#my @preys = ("YY1",	"ZFP41", "ZNF16","ZNF45","ZNF99","ZNF143","ZNF547","ZNF660", "ZNF768","ZNF786",
#			 "ZNF776", "ZNF780A", "ZNF295", "ZNF322", "ZNF397", "ZNF498", "ZNF607", "ZNF112",
#			"ZNF214", "ZNF280C", "PATZ1", "ZFP161"); #"ZNF284", "ZNF708", removed for low expression level

my %exclude = ("ZNF284" => 1, "ZNF708" => 1, "ZNF383" => 1);
my %positiveControl = ("Control4" => 1);
my %negativeControl = ("Control1" => 1, "Control2" => 1, "Control3" => 1, "H2O" => 1, "Flag" => 1, "H2O-GFP" => 1, "Flag-GFP" => 1);
my %all_preys = ();			
my $dir = "./SecondProtocol";
opendir(DIR, $dir) or die $!;
while (my $file = readdir(DIR)) {
	if($file =~ m/^(.+)_NT.*/){
		my $prey = $1;
		unless(exists $exclude{$prey}){
			$all_preys{$prey} = 1;
		}
	}
}
close(DIR);

chdir($dir);
my $result_dir = "Results_mar082018_positiveControl";
mkdir $result_dir;

my @preys = sort(keys %all_preys);
my @treatments = ("NT", "RN_LYSIS", "2NDREAD"); #,


foreach my $treat (@treatments){
	my @existing_preys = ();
	my %prey_bait_avg = ();
	my %prey_bait_stdev = ();
	my %baits = ();
	my $out_treat = $treat;
	$out_treat =~ s/2/S/g;
	my $r_file = $result_dir."/".$out_treat."_replicates_and_average_log_normalized_lum.tab";
	my $r1_file = $result_dir."/".$out_treat."_replicates_only_log_normalized_lum.tab";
	open(R, ">$r_file") or die $!;
	print R "prey\tbait\tmean_lum\tstdev_lum\treplicates\n"; 
	open(R1, ">$r1_file") or die $!;
	print R1 "prey\trep_id\tbait\tlog_lum_zscore\n";
	foreach my $prey (sort @preys){
		my $input = $prey."_".$treat.".tab";
		unless(-e $input){
			print "file $input does not exist\n";
			next;
		}
		push(@existing_preys, $prey);
		my %p384_wellid_gene = (); #p384_wellid as key, gene as value
		my %gene_normalized_lum = (); #gene as key, array of normalized lums as value
		my %gene_AVGlum = (); #gene as key, average of lums as value
		my %p384_wellid_p96_lum = ();	#p384_wellid as key, lum as value, one key per well 
		my %p384_wellid_p96_lum_normalized = ();	#p384_wellid as key, lum as value, one key per well 
		my %p384_p96_positive = (); #p384_p96 as key, array of lums in positive controls as value, one key per 96-well plate
		my %p384_p96_negative = (); #p384_p96 as key, array of lums in negative controls as value, one key per 96-well plate
		my %p384_p96_stats = (); #p384_p96 as key, mean and stdev of negative control lums as value, , one key per 96-well plate
		
		#print "processing $input\n";	
		open(IN, "<$input") or die "can not open $input";
		<IN>;
		<IN>;
		my $blank_count = 0;
		
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
			
			$lum = log($lum);
			if(exists $positiveControl{$gene}){
				push(@{$p384_p96_positive{"$plate384\t$plate96"}}, $lum);
			}elsif(exists $negativeControl{$gene}){
				push(@{$p384_p96_negative{"$plate384\t$plate96"}}, $lum);
			}elsif($gene ne "Blank"){
				$p384_wellid_gene{"$plate384\t$wellID"} = $gene;
				$p384_wellid_p96_lum{"$plate384\t$wellID\t$plate96"} = $lum;
				$baits{$gene} = 1;
			}else{
				$blank_count++;
			}

		}
		print "$treat\t$prey\t$blank_count\n";
		close(IN);
		
		while(my ($key, $nref) = each %p384_p96_negative){
			my @nar = @$nref;
			my @nar_norm = ();
			my $pref = $p384_p96_positive{$key};
			my @par = @$pref;
			
			my $pstat = Statistics::Descriptive::Full->new();
			$pstat->add_data(@par);
			my $pmean = $pstat->mean();
			foreach my $nlum(@nar){
				$nlum /= $pmean;
				push(@nar_norm, $nlum);
			}
			my $nstat = Statistics::Descriptive::Full->new();
			$nstat->add_data(@nar_norm);
			my $nmean = $nstat->mean();
			my $stdev = $nstat->standard_deviation();
			
			my @stats = ($pmean, $nmean, $stdev);
			$p384_p96_stats{$key} = \@stats;
			
			#print "$input\t$key\t$pmean\t$nmean\t$stdev\n";
		}
		
		my $ind = 0;
		while(my ($p384_wellid_p96, $lum) = each %p384_wellid_p96_lum){
			
			my @ids = split(/\t/, $p384_wellid_p96);
			my $p384 = $ids[0];
			my $wellid = $ids[1];
			my $p96 = $ids[2];
			
			my $stats_ref = $p384_p96_stats{"$p384\t$p96"};
			my @stats = @$stats_ref;
			my $pmean = $stats[0];
			my $nmean = $stats[1];
			my $stdev = $stats[2];
			
			unless(defined $pmean){next;}
			my $normalized_lum = 0;
			if($pmean ne 0){
				$lum /= $pmean;
				$normalized_lum = ($lum-$nmean)/$stdev;
			}

			$p384_wellid_p96_lum_normalized{$p384_wellid_p96} = $normalized_lum;
			my $gene = $p384_wellid_gene{"$p384\t$wellid"};
				
			push(@{$gene_normalized_lum{$gene}}, $normalized_lum);
			
			
			$ind++;
			#print "$ind\t$treat\t$prey\t$p384\t$wellid\t$p96\t$gene\t$lum\t$median\t$normalized_lum\n";
			
		}
		
		foreach my $gene (sort(keys %gene_normalized_lum)){
			my $aryref = $gene_normalized_lum{$gene};
			my @l = @$aryref;
			my $stat = Statistics::Descriptive::Full->new();
			$stat->add_data(@l);
			my $mean_lum = $stat->mean();
			#print "before $mean_lum\n";
			$mean_lum =~ s/,//g;
			#print "after $mean_lum\n";
			my $stdev_lum = $stat->standard_deviation();
			$prey_bait_avg{$prey}->{$gene} = $mean_lum;
			$prey_bait_stdev{$prey}->{$gene} = $stdev_lum;
			
			my $sl = join(", ", @l);
			print R "$prey\t$gene\t$mean_lum\t$stdev_lum\t$sl\n";
			for(my $i=0; $i<scalar(@l); $i++){
				my $rep = "rep".$i;
				my $zscore = $l[$i];
				print R1 "$prey\t$rep\t$gene\t$zscore\n";
				
			}
			
		}
		
	}
	
	close(R);
	close(R1);
	
	my $out = $result_dir."/".$out_treat."_bait_prey_log_matrix.tab";
	my $out1 = $result_dir."/".$out_treat."_bait_prey_log_matrix_cluster.tab";
	my @all_baits = sort(keys %baits);
	
	open(OUT, ">$out");
	open(OUT1, ">$out1");
	my $p = join("\tSD\t", @existing_preys);
	my $p1 = join("\t", @existing_preys);
	print OUT "\t$p\n";
	print OUT1 "\t$p1\n";
	foreach my $bait(@all_baits){
		print OUT "$bait";
		print OUT1 "$bait";
		
		foreach my $prey (@existing_preys){
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