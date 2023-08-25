#!/usr/bin/perl -w

#Author: Shuye Pu
#Created on: Feb 10, 2017
#last modified on:

#working directory C:\RSYNC\LUMIER_C2H2\SecondProtocol\Results_mar082018_positiveControl\replicates

#my $lumier = "NT_RNASE_joint_thresholded_network.tab";


my $venn_command2 = "C:/cygwin64/bin/perl C:/RSYNC/worksapce2015/perl_utility/venn_diagram_twoNets.pl";
my $venn_command3 = "C:/cygwin64/bin/perl C:/RSYNC/worksapce2015/perl_utility/venn_diagram_threeNets.pl";

foreach my $dir("Results_mar082018", "Results_mar082018_positiveControl"){
	foreach my $method("replicates", "matrix"){
		my $wd = "C:/RSYNC/LUMIER_C2H2/SecondProtocol/".$dir."/".$method;
		my @nets = ();
		foreach my $treat("NT", "RN_LYSIS", "SNDREAD"){
			my $lumier = $wd."/".$treat."_thresholded_network_".$method.".tab";
			push(@nets, $lumier);
		}
			
		`$venn_command3 $nets[0] $nets[1] $nets[2]`;
		print "done\t$dir\t$method\n";
		
	}
	
	foreach my $treat("NT", "RN_LYSIS", "SNDREAD"){
		my @nets = ();
		foreach my $method("replicates", "matrix"){
			my $wd = "C:/RSYNC/LUMIER_C2H2/SecondProtocol/".$dir."/".$method;
			my $lumier = $wd."/".$treat."_thresholded_network_".$method.".tab";
			push(@nets, $lumier);
		}
			
		`$venn_command2 $nets[0] $nets[1]`;
		print "done\t$dir\t$treat\n";
		
	}
}
foreach my $treat("NT", "RN_LYSIS", "SNDREAD"){
	
	foreach my $method("replicates", "matrix"){
		my @nets = ();
		foreach my $dir("Results_mar082018", "Results_mar082018_positiveControl"){
			my $wd = "C:/RSYNC/LUMIER_C2H2/SecondProtocol/".$dir."/".$method;
			my $lumier = $wd."/".$treat."_thresholded_network_".$method.".tab";
			push(@nets, $lumier);
		}
		`$venn_command2 $nets[0] $nets[1]`;
		print "done\t$treat\t$method\n";
	}
		
}
print "done";