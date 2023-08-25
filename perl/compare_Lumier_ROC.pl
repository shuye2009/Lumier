#!/usr/bin/perl -w

#Author: Shuye Pu
#Created on: Feb 10, 2017
#last modified on:

#working directory C:\RSYNC\LUMIER_C2H2\SecondProtocol\Results_mar082018_positiveControl\replicates

#my $lumier = "NT_RNASE_joint_thresholded_network.tab";


my $ROC_command2 = "C:/cygwin64/bin/perl C:/RSYNC/worksapce2015/MASS_PEC_SCORE/Assessing_PPI_with_ROC.pl";


foreach my $dir("Results_mar082018", "Results_mar082018_positiveControl"){ #
	
	my $wd = "C:/RSYNC/LUMIER_C2H2/SecondProtocol/".$dir."/";
	chdir($wd);
	my @nets = ();
	foreach my $treat("NT", "RN_LYSIS", "SNDREAD"){
		my $lumier = $wd."/".$treat."_replicates_and_average_log_normalized_lum.tab";
		push(@nets, $lumier);
	}
	
	my $inputs = join(",", @nets);
	my $columns = "2,2,2";
	my $outfile = $wd."ROC.png";
	my $labels = "NT,RN_LYSIS,SNDREAD";
		
	`$ROC_command2 $inputs $columns $outfile $labels`;
	print "done\t$dir\n";
	
	
}
print "done";