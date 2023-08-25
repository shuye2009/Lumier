#!/usr/bin/perl -w

#Author: Shuye Pu
#Created on: Feb 10, 2017
#last modified on:

#working directory C:\RSYNC\LUMIER_C2H2\NOV232016
#working directory March 26, 2018 C:\RSYNC\LUMIER_C2H2\SecondProtocol\Results_mar082018_positiveControl

sub max{
	my @input = @{$_[0]};
	my $mx = $input[0];
	for(my $i=1; $i <scalar(@input); $i++){
		if($input[$i] > $mx){
			$mx = $input[$i];
		}
	}
	return $mx;
}

#foreach my $set ("NT", "RNASE"){
foreach my $set ("NT", "RN_LYSIS", "SNDREAD"){	

my $lumir_cutoff = -999;
my $MS_cutoff = 0;
my $saint_cutoff = 0;
my $HG_cutoff = 0;


my $lumier_unthresholded = $set."_replicates_and_average_log_normalized_lum.tab";
my $MS = "C:/RSYNC/LUMIER_C2H2/Nov_10_2016_EM/inter_Nov_10_2016_EM.dat";
my $hgscore = "C:/RSYNC/LUMIER_C2H2/Nov_10_2016_EM/inter_Nov_10_2016_EM_decontaminated_ihgscore0p5.txt";
my $saintscore = "C:/RSYNC/AP_MS_C2H2/PPI_SAINT_2017/Final_data_for_SAINT_analysis/SAINT_score_list_Jan272017.txt";

open(SAINT, "<$saintscore");
open(HG, "<$hgscore");

open(LU, "<$lumier_unthresholded");
open(S, "<$MS");
open (OUTlm, ">Lumier_MS_matching_".$set."_log.tab");
open (OUTml, ">MS_Lumier_matching_".$set."_log.tab");


my %MS_net = ();
my %MS_bait = ();
my %MS_prey = ();
my %common_bait = ();
my %common_prey = ();

while (<S>){
	chomp($_);
    my $temp = $_;
    $temp =~ s/\n//;
    $temp =~ s/\r//;

	
    my @tary = split(/\t/,$temp);
    
    my $baitID = $tary[1];
    my $preyID = $tary[2];
    
    my @baitids = split(/\|/, $baitID);
    my $bait = $baitids[1];
    my @preyids = split(/\|/, $preyID);
    my $prey = $preyids[1];
    
    my $score = $tary[3];
    
    if($score > $MS_cutoff){
    	push(@{$MS_net{"$bait\t$prey"}}, $score);
    	#push(@{$MS_net{"$prey\t$bait"}}, $score);
    	
    	$MS_bait{"$bait"} = 1;
    	$MS_prey{"$prey"} = 1;
    	#print "MS data\t$bait\t$prey\t$score\n";
    }
}

my %hg_net = ();
<HG>;
while (<HG>){
	chomp($_);
    my $temp = $_;
    $temp =~ s/\n//;
    $temp =~ s/\r//;

	
    my @tary = split(/\t/,$temp);
    
    my $baitID = $tary[0];
    my $preyID = $tary[1];
    
    my @baitids = split(/\|/, $baitID);
    my $bait = $baitids[1];
    my @preyids = split(/\|/, $preyID);
    my $prey = $preyids[1];
    
    my $score = $tary[2];
    
    if(defined $prey and defined $bait and $score > $HG_cutoff){
    	
    	$score = sprintf("%.3f", $score);
    	$hg_net{"$bait\t$prey"} = $score;	
    	$hg_net{"$prey\t$bait"} = $score;
    	print "$prey\t$bait\t$score\n";
    }
}

my %saint_net = ();
<SAINT>;
while (<SAINT>){
	chomp($_);
    my $temp = $_;
    $temp =~ s/\n//;
    $temp =~ s/\r//;

	
    my @tary = split(/\t/,$temp);
    
    my $bait = $tary[0];
    my $prey = $tary[1];
    
    my $score = $tary[12];
    
    if($score > $saint_cutoff){
    	$saint_net{"$bait\t$prey"} = $score;
    	$saint_net{"$prey\t$bait"} = $score;
    	print "SAINT\t$prey\t$bait\t$score\n";
    }
}

my %LU_net = ();
my %LU_bait = ();
my %LU_prey = ();
print OUTlm "bait\tprey\tlum_score\tMS_baitPrey\tMS_preyBait\tbaitInMS\tpreyInMS\tsaint\thg\n";
<LU>;
while (<LU>){
	chomp($_);
    my $temp = $_;
    $temp =~ s/\n//;
    $temp =~ s/\r//;

	
    my @tary = split(/\t/,$temp);
    
    my $lm_bait = $tary[0];
    my $lm_prey = $tary[1];
    my $lum_score = $tary[2];
    
    my $bait = $lm_bait;
    my $prey = $lm_prey;
    
    $bait =~ s/"//g;
    $prey =~ s/"//g;
    
    if($bait =~ m/^(.+)?_.*/){
    	$bait = $1;
    }
    if($prey =~ m/^(.+)?_.*/){
    	$prey = $1;
    }
    
    if($prey =~ m/Control|GFP|Blank/g){
    	#print "exclude $prey\n";
    	next;
    }
    my @MS_score1 = ();
    my @MS_score2 = ();
    my $ms_bait = "No";
    my $ms_prey = "No";
    if(exists $MS_net{"$bait\t$prey"}){
    	@MS_score1 = @{$MS_net{"$bait\t$prey"}};
    }
    if(exists $MS_net{"$prey\t$bait"}){
    	@MS_score2 = @{$MS_net{"$prey\t$bait"}};
    }
    if(exists $MS_bait{$bait}){
    	$ms_bait = "MSbait";
    	$common_bait{$bait} = 1;
    }
    if(exists $MS_prey{$bait}){
    	$ms_bait = $ms_bait."MSprey";
    	
    }
    
    if(exists $MS_bait{$prey}){
    	$ms_prey = "MSbait";
    	
    }
    if(exists $MS_prey{$prey}){
    	$ms_prey = $ms_prey."MSprey";
    	$common_prey{$prey} = 1;
    }
    if(exists $common_bait{$bait} and exists $common_prey{$prey}){
    	my $hg = 0;
    	if(exists $hg_net{"$bait\t$prey"}){
    		$hg = $hg_net{"$bait\t$prey"};
    	}
    	
    	my $saint = 0;
    	if(exists $saint_net{"$bait\t$prey"}){
    		$saint = $saint_net{"$bait\t$prey"};
    	}
    	
    	$lum_score = sprintf("%.3f", $lum_score);
    	$hg = sprintf("%.3f", $hg);
    	print OUTlm "$lm_bait\t$lm_prey\t$lum_score\t@MS_score1\t@MS_score2\t$ms_bait\t$ms_prey\t$saint\t$hg\n";
    }
    if($lum_score > $lumir_cutoff){
    	push(@{$LU_net{"$bait\t$prey"}}, $lum_score);
    }
    
    
}

print OUTml "ms_bait\tms_prey\tms_scores\tlum1_baitPrey\tmax1\tlum2_preyBait\tmax2\tmax\tsaint\thg\n";
while(my ($int, $ms_score)= each %MS_net){
	my @ints = split(/\t/, $int);
	my $ms_bait = $ints[0];
	my $ms_prey = $ints[1];
	my @ms_scores = @$ms_score;
	
	my @lum_score1 = ();
	my @lum_score2 = ();
	
	if(exists $LU_net{"$ms_bait\t$ms_prey"}){
		@lum_score1 = @{$LU_net{"$ms_bait\t$ms_prey"}};
	}
	if(exists $LU_net{"$ms_prey\t$ms_bait"}){
		@lum_score2 = @{$LU_net{"$ms_prey\t$ms_bait"}};
	}
	my $max = -999.0;
	my $max1 = max(\@lum_score1);
	if(defined $max1){
		$max = $max1;
	}else{
		$max1 = "";
	}
	
	my $max2 = max(\@lum_score2);
	if(defined $max2){
		if($max2 > $max){
			$max = $max2;
		}
	}else{
		$max2 = "";
	}
	
	if(exists $common_bait{$ms_bait} and exists $common_prey{$ms_prey}){
		my $hg = 0;
    	if(exists $hg_net{"$ms_bait\t$ms_prey"}){
    		$hg = $hg_net{"$ms_bait\t$ms_prey"};
    	}
    	
    	my $saint = 0;
    	if(exists $saint_net{"$ms_bait\t$ms_prey"}){
    		$saint = $saint_net{"$ms_bait\t$ms_prey"};
    	}
    	print OUTml "$ms_bait\t$ms_prey\t@ms_scores\t@lum_score1\t$max1\t@lum_score2\t$max2\t$max\t$saint\t$hg\n";
    }
}

my @cb = keys %common_bait;

my $ncb = scalar(@cb);

#print OUT "\nNumber of common baits\t$ncb\n";

}
print "done";