#!/usr/bin/perl -w

#Author: Shuye Pu
#Created on: Oct 18, 2016
#last modified on: March 08, 2018

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
my $dir = "./Archive1_3";
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
my $result_dir = "Results_oct032018_same_median_normalization";
mkdir $result_dir unless (-e $result_dir and -d $result_dir);

my @preys = sort(keys %all_preys);
my @treatments = ("NT", "RN_LYSIS", "2NDREAD"); #,
my %out_treatments = ("NT" => "NT", "RN_LYSIS" => "RN", "2NDREAD" => "SNDREAD");
my $max_lum = 7488349; ## determed by preprocessing
my %prey_treat_384_96_median = (); # store median for each 96-well plate to get median of median, which will be used to scale each 96-well plate
my @prey_treat_384_96_median_array = ();
my $median_file = $result_dir."/median_of_all_96well_plates.tab";
open(M, ">$median_file") or die $!;
foreach my $treat (@treatments){
	my @existing_preys = ();
	my %prey_bait_avg = ();
	my %prey_bait_stdev = ();
	my %baits = ();
	my $out_treat = $out_treatments{$treat};
	
	
	foreach my $prey (sort @preys){
		my $input = $prey."_".$treat.".tab";
		unless(-e $input){
			print "file $input does not exist\n";
			next;
		}
		push(@existing_preys, $prey);
		my %p384_p96_lum = (); #p384_p96 as key, array of lums as value, one key per 96-well plate
		
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
			
			
			if($lum eq "OVRFLW"){
				$lum = $max_lum;
			}elsif($lum eq "X"){
				$lum = -1;
			}
			if($lum > $max_lum){
				$max_lum = $lum;
			}
			
			### only consider valid lum values
			if(defined $gene and $lum > 0){
				push(@{$p384_p96_lum{"$plate384\t$plate96"}}, $lum);
			}	
			
		}
		
		close(IN);
		
		while(my ($key, $ref) = each %p384_p96_lum){
			my @ar = @$ref;
			my $stat = Statistics::Descriptive::Full->new();
			$stat->add_data(@ar);
			my $median = $stat->median();
			push(@prey_treat_384_96_median_array, $median);
			$prey_treat_384_96_median{$prey}{$treat}{$key} = $median;
			print M "$treat\t$prey\t$key\t$median\n";
		}
	}
}
close(M);
my $stat = Statistics::Descriptive::Full->new();
$stat->add_data(@prey_treat_384_96_median_array);
my $median_of_median = $stat->median();

my %constructs_length = ();
my %prey_treat_384_well_gene_lum = ();
my %prey_treat_384_well_gene_nlum = (); # store all values to output tables for DESeq and HGscore
foreach my $treat (@treatments){
	my @existing_preys = ();
	my %prey_bait_avg = ();
	my %prey_bait_stdev = ();
	my %baits = ();
	my $out_treat = $out_treatments{$treat};
	my $r_file = $result_dir."/".$out_treat."_replicates_and_average_median_normalized_lum.tab";
	my $r1_file = $result_dir."/".$out_treat."_replicates_only_median_normalized_lum.tab";
	open(R, ">$r_file") or die $!;
	print R "prey\tbait\tmean_lum\tstdev_lum\treplicates\n"; 
	open(R1, ">$r1_file") or die $!;
	print R1 "prey\trep_id\tbait\tmedian_normalized_lum\n";
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
					$constructs_length{$gene} = 300;
				}
				
			}
			
			
			#unless(defined $gene){$gene = "Blank";}
			
			if($lum eq "OVRFLW"){
				$lum = $max_lum;
			}elsif($lum eq "X"){
				$lum = -1;
			}
			
			if(defined $gene){
				$baits{$gene} = 1;
				$all_baits{$gene} = 1;
			}else{
				$blank_count++;
				$gene = "Blank";
			}
			if($lum > 0){
				$p384_wellid_gene{"$plate384\t$wellID"} = $gene;
				$p384_wellid_p96_lum{"$plate384\t$wellID\t$plate96"} = $lum;
			}
		}
		
		close(IN);
		
		
		my $ind = 0;
		while(my ($p384_wellid_p96, $lum) = each %p384_wellid_p96_lum){
			
			my @ids = split(/\t/, $p384_wellid_p96);
			my $p384 = $ids[0];
			my $wellid = $ids[1];
			my $p96 = $ids[2];
			
			my $median = $prey_treat_384_96_median{$prey}{$treat}{"$p384\t$p96"};
			
			unless(defined $median){next;}
			my $scale_factor = $median/$median_of_median;
			my $normalized_lum = 0;
			if($scale_factor ne 0){
				$normalized_lum = $lum/$scale_factor;
			}

			#print "$input\t$p384_wellid_p96\t$lum\t$q1\t$median\t$q3\t$stdev\t$normalized_lum\n";
			$p384_wellid_p96_lum_normalized{$p384_wellid_p96} = $normalized_lum;
			my $gene = $p384_wellid_gene{"$p384\t$wellid"};
				
			push(@{$gene_normalized_lum{$gene}}, $normalized_lum);
			
			$prey_treat_384_well_gene_lum{$gene}{$prey}{$treat}{"$p384\t$p96\t$wellid"} = $lum;
			$prey_treat_384_well_gene_nlum{$gene}{$prey}{$treat}{"$p384\t$p96\t$wellid"} = $normalized_lum;
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
				my $nlum = $l[$i];
				print R1 "$prey\t$rep\t$gene\t$nlum\n";
				
			}
			
		}
		print "$treat\t$prey\t$blank_count\n";
	}
	
	close(R);
	close(R1);
	
	my $out = $result_dir."/".$out_treat."_bait_prey_matrix.tab";
	my $out1 = $result_dir."/".$out_treat."_bait_prey_matrix_cluster.tab";
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

my $OR_out = $result_dir."/table_for_original_lum.tab";
my $DE_out = $result_dir."/table_for_differential_analysis.tab";
my $COMPASS_out_NT = $result_dir."/NT_table_for_COMPASSscore.tab";
my $COMPASS_out_RNASE = $result_dir."/RN_table_for_COMPASSscore.tab";
my $COMPASS_out_SNDREAD = $result_dir."/SNDREAD_table_for_COMPASSscore.tab";

open(OR, ">$OR_out") or die $!; # original lum values
open(DE, ">$DE_out") or die $!; # median normalized lum values
open(RNASE, ">$COMPASS_out_RNASE") or die $!;
open(NT, ">$COMPASS_out_NT") or die $!;
open(SNDREAD, ">$COMPASS_out_SNDREAD") or die $!;
foreach my $fh (NT, RNASE, SNDREAD){
	print $fh "Purification\tBait\tPrey\tSpectral Counts\n";
}
my %fh_hash = ("NT" => NT, "RN_LYSIS" => RNASE, "2NDREAD" => SNDREAD);

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
			my $fh = $fh_hash{$treat};
			my $rep1 = "NA";
			my $rep2 = "NA";
			my $rep1or = "NA";
			my $rep2or = "NA";
			my $loc1 = "NA";
			my $loc2 = "NA";
			if(defined $loc_hash_ref){
				my %loc_hash = %$loc_hash_ref;
				my %loc_hash_or = %$loc_hash_ref_or;
				my @locs = sort keys %loc_hash;
				
				if(defined $locs[0]){
					$loc1 = $locs[0];
					$rep1 = $loc_hash{$locs[0]};
					$rep1or = $loc_hash_or{$locs[0]};
				}
				if(defined $locs[1]){
					$loc2 = $locs[1];
					$rep2 = $loc_hash{$locs[1]};
					$rep2or = $loc_hash_or{$locs[1]};
				}
				## output for COMAPSS score
				foreach my $loc (@locs){
					my $rep = $loc_hash{$loc};
					$loc =~ s/\t/_/g;
					if(defined $rep) {print $fh "$loc\t$bait\t$prey\t$rep\n";} #$rep = log($rep)/log(10);
				}
				
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
print "finish!\nmedian of medians: $median_of_median\nMax lum: $max_lum\n";