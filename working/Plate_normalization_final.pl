#!/usr/bin/perl -w

#Author: Shuye Pu
#Created on: Oct 18, 2016
#last modified on: March 08, 2018
# Original file name "Plate_normalization_by_control_distribution.pl"
# This is the final version for Lumier score normalization,
# the normalized zscore is used for deriving interactions

#working directory C:\RSYNC\LUMIER_C2H2

use Statistics::Descriptive qw(:all);
use Regexp::Common;

my $pl_file = "C:/RSYNC/AP_MS_C2H2/human_name_length.txt";
open(PL, "<$pl_file") or die $!;
my %length_hash = ();
while(<PL>){
	chomp($_);
	my $temp = $_;
    $temp =~ s/\n//;
    $temp =~ s/\r//;
    my @tary = split(/\t/,$temp);
    $length_hash{$tary[0]} = $tary[1];
}

my %exclude = ();
#my %exclude = ("ZNF284" => 1, "ZNF708" => 1, "ZNF383" => 1);
my %all_preys = ();	
my %all_baits = ();			
#my $dir = "./SecondProtocol";
my $dir = "./Archive1_3"; # "./Archive1_3";
opendir(DIR, $dir) or die $!;
while (my $file = readdir(DIR)) {
	if($file =~ m/^(.+)_NT.*\.tab/){
		my $prey = $1;
		unless(exists $exclude{$prey}){
			$all_preys{$prey} = 1;
		}
	}
}
close(DIR);

chdir($dir);
my $result_dir = "Results_oct032018_control_distribution_normalization";
mkdir $result_dir unless (-e $result_dir and -d $result_dir);

my %prey_controls = ("Control1" => 1,"Control2" => 1,"Control3" => 1,"Control4" => 1); ### seems not used
my %bait_controls = ("Flag" => 1, "Flag-GFP" => 1, "H2O" => 1, "H2O-GFP" => 1);  # Blank is not cosidered as it is missing in same plate

my @preys = sort(keys %all_preys);
my @treatments = ("NT", "RN_LYSIS", "2NDREAD"); #,"DNAse"); # 
my %out_treatments = ("NT" => "NT", "RN_LYSIS" => "RN", "2NDREAD" => "SNDREAD"); #"DNAse" => "DNASE"); #
my $max_lum = 7441931; ## determined by preprocessing
my %treat_baitControl_mean = (); # store mean of prey controls for each treat, which will be used to scale each 96-well plate
my %treat_baitControl_stdev = ();
my %treat_prey_plate_mean = (); # store median of each 96plate for each treat and each prey, which will be used to scale each 96-well plate


## collect control stats for bait and prey
foreach my $treat (@treatments){
	my $treato = $out_treatments{$treat};
	my $control_file = $result_dir."/".$treato."_control_of_all_96well_plates.tab";
	open(M, ">$control_file") or die $!;

	## get means for prey controls and bait controls
	print M "Median\tMean\tSTDEV\tMIN\tMAX\tFlag\tH2O\n";
	
	my %prey_control_lum = (); # prey_control as key, array of lums as value 
	my %bait_control_lum = (); #bait control as key, array of lums as value
	foreach my $prey (sort @preys){
		my $input = $prey."_".$treat.".tab";
		unless(-e $input){
			print "file $input does not exist\n";
			next;
		}
		
		my %p384_p96_control_lum = (); #p384_p96 and control type as key, array of lums as value, one key per 96-well plate
		
		
		my %p384_p96_lum = (); #p384_p96, array of lums as value, one key per 96-well plate
		#print "processing $input\n";	
		open(IN, "<$input") or die "can not open $input";
		
		while(<IN>){
			chomp($_);
		    my $temp = $_;
		    $temp =~ s/\n//;
		    $temp =~ s/\r//;
					
		    my @tary = split(/\t/,$temp);
		    
		    #if($input eq "JARID1B_RNASE.tab"){print "$input\t$size\t$temp\n";}
		    unless(defined $tary[0]){next;}
			if($tary[0] eq "" or $tary[0] eq "Plate Number"){next;}
			
			my $plate384 = $tary[0];
			
			my $lum = $tary[3];
			my $plate96 = $tary[4];
			my $gene = $tary[7];
			
			
			unless(defined $gene){
				$gene = "Blank";
			}
			if($gene eq "Flag-GFP"){
				$gene = "Flag";
			}elsif($gene eq "H2O-GFP"){
				$gene = "H2O";
			}
			
			if($lum eq "OVRFLW"){
				$lum = $max_lum;
			}elsif($lum eq "X"){
				$lum = 1;
			}
			if($lum > $max_lum){
				$max_lum = $lum;
			}
			$lum = log($lum);
			### only consider valid lum values
			if($lum > 0){
				if(exists $prey_controls{$gene}){
					push(@{$p384_p96_control_lum{"$plate384\t$plate96"}->{$gene}}, $lum);
					push(@{$prey_control_lum{$gene}}, $lum);
				}elsif(exists $bait_controls{$gene}){
					push(@{$p384_p96_control_lum{"$plate384\t$plate96"}->{$gene}}, $lum);
					push(@{$bait_control_lum{$gene}}, $lum);
				}else{
					push(@{$p384_p96_lum{"$plate384\t$plate96"}}, $lum);
				}
				
			}	
			
		}
		
		close(IN);
		
		# output all control lum values for inspection purpose
		
		while(my ($key, $href) = each %p384_p96_control_lum){
			my @averages_controls = ();
			my @typeof_controls = ();
			my %control_hash = %$href;
			foreach my $c (sort keys %control_hash){
				my $ref = $control_hash{$c};
				my @ar = @$ref;
				my $stat = Statistics::Descriptive::Full->new();
				$stat->add_data(@ar);
				my $mean = $stat->mean();
				push(@averages_controls, $mean);
				push(@typeof_controls, $c);
			}
			my $array_ref = $p384_p96_lum{$key};
			my @ar = @$array_ref;
			my $stat = Statistics::Descriptive::Full->new();
			$stat->add_data(@ar);
			my $median = $stat->median();
			my $mean = $stat->mean();
			my $std = $stat->standard_deviation();
			my $min = $stat->min();
			my $max = $stat->max();
			$treat_prey_plate_mean{$treat}{$prey}{$key} = $mean;
			
			my $individual_controls = join("\t", @averages_controls);
			my $row_id = "$treat\t$prey\t$key";
			$row_id =~ s/\t/_/g;
			$row_id =~ s/\s//g;
			print M "$row_id\t$median\t$mean\t$std\t$min\t$max\t$individual_controls\n";
		}
		
	}
	while(my ($key, $ref) = each %bait_control_lum){
		my @ar = @$ref;
		my $stat = Statistics::Descriptive::Full->new();
		$stat->add_data(@ar);
		my $mean = $stat->mean();
		my $std = $stat->standard_deviation();
		$treat_baitControl_mean{$treat}{$key} = $mean;
		$treat_baitControl_stdev{$treat}{$key} = $std;
		#print M "$treat\t$key\t$mean\n";
	}
	close(M);
}




## output normalized lum for downstream analysis
if(1){
my %constructs_length = ();
my %prey_treat_384_well_gene_lum = ();
my %prey_treat_384_well_gene_nlum = (); # store all values to output tables for DESeq and HGscore
foreach my $treat (@treatments){
	
	my $control = "Flag";
	my $control_mean = $treat_baitControl_mean{$treat}{$control};
	my $control_stdev = $treat_baitControl_stdev{$treat}{$control};
	
	my @existing_preys = ();
	my %prey_bait_avg = ();
	my %prey_bait_max = ();  # output for ROC evaluation
	my %prey_bait_stdev = ();
	my %baits = ();
	my $out_treat = $out_treatments{$treat};
	my $r_file = $result_dir."/".$out_treat."_replicates_and_average_control_normalized_lum.tab";
	my $r1_file = $result_dir."/".$out_treat."_replicates_only_control_normalized_lum.tab";
	open(R, ">$r_file") or die $!;
	print R "prey\tbait\tmean_lum\tstdev_lum\treplicates\n"; 
	open(R1, ">$r1_file") or die $!;
	print R1 "prey\trep_id\tbait\tcontrol_normalized_lum\n";
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
		my %p384_preyControl_p96_lum = ();	#p384_p96_preyControl as key, lum as value, one key per control
		my %p384_baitControl_p96_lum = ();	#p384_p96_baitControl as key, lum as value, one key per control
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
				my $len = $length_hash{$gene};
				$gene = $gene."_".$construct;
				if(defined $len){
					$constructs_length{$gene} = $len;
				}else{
					$constructs_length{$gene} = 300; # average protein length
				}
				
			}
			
			
			
			if(defined $gene){
				if($gene eq "Flag"){
					$gene = "Flag-GFP";
				}elsif($gene eq "H2O"){
					$gene = "H2O-GFP";
				}
				$baits{$gene} = 1;
				$all_baits{$gene} = 1;
			}else{
				$blank_count++;
				$gene = "Blank";
			}
			
			
			if($lum eq "OVRFLW"){
				$lum = $max_lum;
			}elsif($lum eq "X"){
				$lum = 1;
			}
			# $lum = log($lum);
			
			
			if($lum > 0){
				$p384_wellid_gene{"$plate384\t$plate96"}{$wellID} = $gene;
				$p384_wellid_p96_lum{"$plate384\t$plate96"}{$wellID} = $lum;
				if(exists $prey_controls{$gene}){
					$p384_preyControl_p96_lum{"$plate384\t$plate96"}{$gene} = $lum;
				}elsif(exists $bait_controls{$gene}){
					$p384_baitControl_p96_lum{"$plate384\t$plate96"}{$gene} = $lum;
				}
			}
		}
		
		close(IN);
		
		
		my $ind = 0;
		while(my ($p384_p96, $lum_href) = each %p384_wellid_p96_lum){
		
			my %lum_hash = %$lum_href;
			while (my ($wellID, $lum) = each %lum_hash){

				my $normalized_lum = 0;
				if($control_stdev ne 0){
					$normalized_lum = (log($lum) - $control_mean)/$control_stdev;
					 
				}
	
				#print "$input\t$p384_wellid_p96\t$lum\t$q1\t$median\t$q3\t$stdev\t$normalized_lum\n";
				$p384_wellid_p96_lum_normalized{"$p384_p96\t$wellID"} = $normalized_lum;
				my $gene = $p384_wellid_gene{$p384_p96}{$wellID};
					
				push(@{$gene_normalized_lum{$gene}}, $normalized_lum);
				
				$prey_treat_384_well_gene_lum{$gene}{$prey}{$treat}{"$p384_p96\t$wellID"} = $lum;
				$prey_treat_384_well_gene_nlum{$gene}{$prey}{$treat}{"$p384_p96\t$wellID"} = $normalized_lum;
				$ind++;
				#print "$ind\t$treat\t$prey\t$p384\t$wellid\t$p96\t$gene\t$lum\t$median\t$normalized_lum\n";
			}
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
			
			my @geneConstruct = split(/_/,$gene);
			my $geneName = $geneConstruct[0];
			if(exists $prey_bait_max{$prey}->{$geneName}){
				my $max = $prey_bait_max{$prey}->{$geneName};
				if ($mean_lum > $max){
					$prey_bait_max{$prey}->{$geneName} = $mean_lum;
				}
			}else{
				$prey_bait_max{$prey}->{$geneName} = $mean_lum;
			}
			
			my $sl = join(", ", @l);
			print R "$prey\t$gene\t$mean_lum\t$stdev_lum\t$sl\n";
			for(my $i=0; $i<scalar(@l); $i++){
				my $rep = "rep".$i;
				my $nlum = $l[$i];
				print R1 "$prey\t$rep\t$gene\t$nlum\n";
				
			}
			
		}
		print "$treat\t$prey\t$blank_count\n";
	}
	
	close(R);
	close(R1);
	
	my $outroc = $result_dir."/".$out_treat."_bait_prey_zscore_network.tab";
	my $out = $result_dir."/".$out_treat."_bait_prey_matrix.tab";
	my $out1 = $result_dir."/".$out_treat."_bait_prey_matrix_cluster.tab";
	my @all_baits = sort(keys %baits);
	
	open(OUTroc, ">$outroc");
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
	
	print OUTroc "Bait\tPrey\tZscore\n";
	foreach my $prey (keys %prey_bait_max){
		my $hashref = $prey_bait_max{$prey};
		my %bait_hash = %$hashref;
		while( my ($bait, $lum_max) = each %bait_hash){
			$bait = uc($bait);
			print OUTroc "$bait\t$prey\t$lum_max\n";
		}
	}
}



my $OR_out = $result_dir."/table_for_original_lum.tab";
my $DE_out = $result_dir."/table_for_zscore.tab";


open(OR, ">$OR_out") or die $!; # original lum values
open(DE, ">$DE_out") or die $!; # zscore normalized lum values

#print DE "Bait";
my @header = ();
foreach my $prey(keys %all_preys){
	foreach my $treat (@treatments){
		my $out_treat = $out_treatments{$treat};
		push(@header, $prey."_".$out_treat."_1");
		push(@header, $prey."_".$out_treat."_2");
	}
}
my $head = join("\t", @header);
print DE "$head\n";
print OR "$head\n";

foreach my $bait (keys %all_baits){
	
	print DE "$bait";
	print OR "$bait";
	foreach my $prey (keys %all_preys){
		foreach my $treat (@treatments){
			my $loc_hash_ref = $prey_treat_384_well_gene_nlum{$bait}{$prey}{$treat};
			my $loc_hash_ref_or = $prey_treat_384_well_gene_lum{$bait}{$prey}{$treat};
			
			my $rep1 = 0; # default zscore
			my $rep2 = 0;
			my $rep1or = exp($treat_baitControl_mean{$treat}{"H2O"}); # default lum score
			my $rep2or = exp($treat_baitControl_mean{$treat}{"H2O"});
			
			if(defined $loc_hash_ref){
				my %loc_hash = %$loc_hash_ref;
				my %loc_hash_or = %$loc_hash_ref_or;
				my @locs = sort keys %loc_hash;
				
				if(defined $locs[0]){
					
					$rep1 = $loc_hash{$locs[0]};
					$rep1or = $loc_hash_or{$locs[0]};
				}
				if(defined $locs[1]){
					
					$rep2 = $loc_hash{$locs[1]};
					$rep2or = $loc_hash_or{$locs[1]};
				}
				
			}else{
				print "$bait\t$prey\t$treat lum is not defined\n";
			}
			print DE "\t$rep1\t$rep2";
			print OR "\t$rep1or\t$rep2or";
			
		}
	}
	print DE "\n";
	print OR "\n";
}

while (my ($c, $l) = each %constructs_length){
	print "$c\t$l\n";
}
}
print "finish!\nMax lum: $max_lum\n";