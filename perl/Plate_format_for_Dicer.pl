#!/usr/bin/perl -w

#Author: Shuye Pu
#Created on: Oct 18, 2016
#last modified on: March 08, 2018

#working directory C:\RSYNC\LUMIER_C2H2

use Statistics::Descriptive qw(:all);
use Regexp::Common;



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
my $result_dir = "For_Dicer";
mkdir $result_dir unless (-e $result_dir and -d $result_dir);

my @preys = sort(keys %all_preys);
my @treatments = ("NT", "RN_LYSIS", "2NDREAD"); #,
my %out_treatments = ("NT" => "NT", "RN_LYSIS" => "RN", "2NDREAD" => "SNDREAD");
my $max_lum = 7488349; ## determed by preprocessing

foreach my $treat (@treatments){
	
	foreach my $prey (sort @preys){
		my $input = $prey."_".$treat.".tab";
		unless(-e $input){
			print "file $input does not exist\n";
			next;
		}
		my $p384_plate1 = $result_dir."/".$prey."_".$treat."_3841.csv";
		my $p384_plate2 = $result_dir."/".$prey."_".$treat."_3842.csv";
		my $p384_delete = $result_dir."/".$prey."_".$treat."_384delete.txt";
		
		my %p384_delete = ();
		my %p384_p96_delete = ();
		open(OUT1, ">$p384_plate1") or die $!;
		open(OUT2, ">$p384_plate2") or die $!;
		open(D, ">$p384_delete") or die $!;
		#print "processing $input\n";	
		open(IN, "<$input") or die "can not open $input";
		
		my $first = 0;
		while(<IN>){
			chomp($_);
		    my $temp = $_;
		    $temp =~ s/\n//;
		    $temp =~ s/\r//;
					
		    my @tary = split(/\t/,$temp);
		    
		    #if($input eq "JARID1B_RNASE.tab"){print "$input\t$size\t$temp\n";}
		    unless(defined $tary[0]){next;}
			if($tary[0] eq "" ){next;}
			if($tary[0] eq "Plate Number"){
				my $header = join(",", @tary);
				$first = $first + 1;
				if($first eq 1){print OUT1 "$header\n";}else{print OUT2 "$header\n";}
				next;
			}else{
				my $line = join(",", @tary);
				my $plate384 = $tary[0];
				my $well = $tary[2];
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
			
				unless(defined $gene){
					push(@{$p384_p96_delete{$plate384}{$plate96}}, $well);
					push(@{$p384_delete{$plate384}}, $well);
				}
				if($first eq 1){print OUT1 "$line\n";}else{print OUT2 "$line\n";}
			}
			
		}
		close(IN);
		foreach my $p384 (sort keys %p384_p96_delete){
			my $p96_ref = $p384_p96_delete{$p384};
			my %p96_hash = %$p96_ref;
			foreach my $p96 (sort keys %p96_hash){
				my $array_ref = $p96_hash{$p96};
				my @deletes = sort(@$array_ref);
				my $deleteWells = join(",", @deletes);
				print D "$p384\t$p96\t$deleteWells\n";
			}
			my $array_ref = $p384_delete{$p384};
			my @deletes = sort(@$array_ref);
			my $deleteWells = join(",", @deletes);
			print D "$p384\t$deleteWells\n";
		}
	}
}
		
print "finish!\nMax lum: $max_lum\n";