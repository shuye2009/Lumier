#!/usr/bin/perl -w

#Author: Shuye Pu
#Created on: Feb 10, 2017
#last modified on: Oct 24, 2018
my %scores = ("maxquant" => {"untreated" => 0.67, "Rnase" => 0.74}); #determined in format_prohits_for_COMPASS.R
foreach my $subdir ("maxquant"){
foreach my $dataset ("untreated", "Rnase"){

my $HGscore_cutoff = $scores{$subdir}{$dataset};
my $HGscore = "C:/RSYNC/LUMIER_C2H2/Archive_oct032018/Archive1_3/Results_oct032018_control_distribution_normalization/HGscore_from_APMS_180913/".$subdir."/".$dataset."_ihgscore0p5.txt";

my $venn_command = "C:/cygwin64/bin/perl C:/RSYNC/worksapce2015/perl_utility/venn_diagram_twoNets.pl";

open(S, "<$HGscore");
my %HGscore_net = ();
#my %HGscore_bait = ();
my %HGscore_prey = ();

<S>;
while (<S>){
	chomp($_);
    my $temp = $_;
    $temp =~ s/\n//;
    $temp =~ s/\r//;

	
    my @tary = split(/\t/,$temp);
    
    my $bait = $tary[0];
    my $prey = $tary[1];
    
    $bait =~ s/\s//g;
    $bait =~ s/"//g;
    $bait = uc($bait);
    $prey =~ s/\s//g;
    $prey =~ s/"//g;
    $prey = uc($prey);
    
    my $score = $tary[2];
    
    if($score > $HGscore_cutoff){
    	if($bait lt $prey){
    		$HGscore_net{"$bait\t$prey"} = $score;
    	}else{
    		$HGscore_net{"$prey\t$bait"} = $score;
    	}
    	
    	
    	#$HGscore_bait{"$bait"} = 1;
    	$HGscore_prey{"$prey"} = 1;
    	$HGscore_prey{"$bait"} = 1;
    	#print "HGscore\t$bait\t$prey\t$score\n";
    	
    }
}

my @ppi = keys %HGscore_net;
my $size_net = scalar(@ppi);
#foreach my $dir ("Results_mar082018", "Results_mar082018_positiveControl"){
#	foreach my $method("replicates", "matrix"){
#		my $wd = "C:/RSYNC/LUMIER_C2H2/SecondProtocol/".$dir."/".$method;
		my $outdir = "C:/RSYNC/LUMIER_C2H2/Archive_oct032018/Archive1_3/Results_oct032018_control_distribution_normalization/HGscore_from_APMS_180913/".$subdir;
		
		my $method = "replicates";
		my $wd = "C:/RSYNC/LUMIER_C2H2/Archive_oct032018/Archive1_3/Results_oct032018_control_distribution_normalization/replicates_min";
		#foreach my $treat("NT", "RN", "SNDREAD"){
			my $treat = "NT";
			if($dataset eq "Rnase"){
				$treat = "RN";
			}
			my $lumier = $wd."/".$treat."_thresholded_network_".$method.".tab";
			
			open(L, "<$lumier");
			my $lumier_basedon_HGscore = $outdir."/".$treat."_".$method."_LUMIER_network.tab";
			open (OUT1, ">$lumier_basedon_HGscore");
			my $HGscore_basedon_lumier = $outdir."/".$treat."_".$method."_".$dataset."_HGSCORE_network_cutoff_".$HGscore_cutoff.".tab";
			open (OUT2, ">$HGscore_basedon_lumier");
			
			
			my %lumier_net = ();
			my %lumier_bait = ();
			my %lumier_prey = ();
			
			my $lumierINhg = 0;
			
			my $h = <L>;
			
			while (<L>){
				chomp($_);
			    my $temp = $_;
			    $temp =~ s/\n//;
			    $temp =~ s/\r//;
			
				
			    my @tary = split(/\t/,$temp);
			    
			    my $bait = $tary[1];
			    my $prey = $tary[0];
			    my $score = $tary[2];
			    
			    $bait =~ s/\s//g;
			    $bait =~ s/"//g;
			    $bait = uc($bait);
			    $prey =~ s/\s//g;
			    $prey =~ s/"//g;
			    $prey = uc($prey);
			    
			    if($bait =~ m/^(.+)?_.*/){
			    	$bait = $1;
			    }
			    if($prey =~ m/^(.+)?_.*/){
			    	$prey = $1;
			    }
			    $lumier_bait{$bait} = 1;
			    $lumier_prey{$prey} = 1;
			    $lumier_net{"$bait\t$prey"} = $score;
			    
			    if(exists $HGscore_prey{$bait} and $HGscore_prey{$prey}){
			    	#if($bait ne $prey){
			    		print OUT1 "$bait\t$prey\t$score\n";
			    		$lumierINhg++;
			    	#}
			    }
			}
			my @lumier_ppi = keys %lumier_net;
			my $lumier_size = scalar(@lumier_ppi);
			close(L);
			close(OUT1);
			my $hgINlumier = 0;
			while(my ($int, $score) = each %HGscore_net){
				my @interactors = split(/\t/, $int);
				my $bait = $interactors[0];
				my $prey = $interactors[1];
				
				if((exists $lumier_bait{$bait} and exists $lumier_prey{$prey}) or (exists $lumier_bait{$prey} and exists $lumier_prey{$bait})){
					print OUT2 "$int\t$score\n";
					$hgINlumier++;
				}
			}
			close(OUT2);
			`$venn_command $lumier_basedon_HGscore $HGscore_basedon_lumier`;
			unlink $lumier_basedon_HGscore;
			unlink $HGscore_basedon_lumier;
			print "done\t$method\t$treat\tLUMIER size: $lumier_size\tLUMIER in HG: $lumierINhg\tHG in LUMIER: $hgINlumier\n";
		#}
	#}
#}
print "done\tSize of compared net: $size_net\n";
}
}