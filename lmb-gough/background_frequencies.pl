#!/usr/bin/perl -w
use strict;

#this is to generate the background.vcf.dat file from a background.vcf
#but restricted to an input list of positions, e.g. 23andMe
#this is *not* the script that is part of Snowflake, the Snowflake script only reports protein-coding regiosn whereas this will do the whole genome, however the snowflake script does not report positions that are 100% same as reference

my $reference="/home/gough/data/snowflake/Homo_sapiens.GRCh37.75.dna.toplevel.fa";
my $background="/home/gough/data/snowflake/1000g/";
#my $api="/home/gough/data/snowflake/api.23andme.com";
#my $eoutput="/home/gough/data/snowflake/background_euro_23andme.dat";
#my $output="/home/gough/data/snowflake/background_23andme.dat";
my $api="/home/gough/data/rihab/ancestryDNA/ancestrydna.list";
my $eoutput="/home/gough/data/snowflake/background_euro_ancestrydna.dat";
my $output="/home/gough/data/snowflake/background_ancestrydna.dat";
my $euro="/home/gough/data/snowflake/igsr_samples.tsv";
my (%euro,%alleles2,%positions1,%positions2,%euno,%outno);

#load in the 23andme positions
open F,("$api");
while (<F>){
  if (/\S/){
    unless (/^#/ or /^index/){
      my @columns=split  /\t/,$_;
      chomp($columns[3]);
      my $id="$columns[2]:$columns[3]";
      $positions1{$id}=1;
    }
  }
}
close F;
print STDERR "loaded ",scalar(keys(%positions1))," 23andme positions\n";


#load in list of Europeans
open E,("$euro");
while (<E>){
  if (/^([A-Z]+\d+)\s/){
    $euro{$1}=1;
  }
  else{
    unless (/^Sample/){
      print STDERR $_;
    }
  }
}
close E;

#load in the background allele frequencies
open EU,(">$eoutput");
open OUT,(">$output");
my @list=("11","10","00","01","1","0");
open LS,("ls  $background |");
while (<LS>){
  if (/^ALL.chr(\w+)\./){
    my $chr=$1;my %e;
    open F,("gunzip -c $background/$_ |");
    while (<F>){
      if (/^#CHROM/){
	my @temp=split /\t/,$_;
	$euno{$chr}=scalar(keys(%euro));$outno{$chr}=scalar(@temp)-9;
	foreach my $i (0 .. scalar(@temp)-1){
	  if (exists($euro{$temp[$i]})){
	    $e{$i}=$temp[$i];
	  }
	}
      }
      else{
	unless (/^#/){
	  my @temp=split /\t/,$_;
	  chomp($temp[scalar(@temp)-1]);
	  if (exists($positions1{"$temp[0]:$temp[1]"})){
	    $positions1{"$temp[0]:$temp[1]"}=0;
	    if ($temp[4] =~ /^([ATCG])\,/){$temp[4]=$1;}
	    if ($temp[4] =~ /^[ATCG]$/ and $temp[3] =~ /^[ATCG]$/){
	      
	      my %eu;my $x=0;
	      foreach my $i (keys(%e)){
		if ($temp[$i] =~ /([0-9\.])\S+([0-9\.])/){
		  $eu{"$1$2"}++;
		}
		elsif($temp[$i] =~ /([0-9\.])/){
		  $eu{"$1"}++;
		}
		else{
		  print STDERR "ERROR: $temp[$i] $temp[0]:$temp[1]\n";
		}
	      }
	      foreach my $a (@list){
		if (exists($eu{$a})){
		  if ($x == 1){
		    print EU ";";
		  }
		  else{
		    print EU "$temp[0]_$temp[1]_$temp[3]/$temp[4]\t";
		  }
		  print EU "$a=$eu{$a}";
		  $x=1;
		}
	      }
	      if ($x == 1){print EU "\n";}
	      
	      
	      my %out;$x=0;
	      foreach my $i (9 .. scalar(@temp)-1){
		if ($temp[$i] =~ /([0-9\.])\S+([0-9\.])/){
		  $out{"$1$2"}++;
		}
		elsif($temp[$i] =~ /([0-9\.])/){
		  $out{"$1"}++;
		}
		else{
		  print STDERR "ERROR: $temp[$i] $temp[0]:$temp[1]\n";
		}
	      }
	      foreach my $a (@list){
		if (exists($out{$a})){
		  if ($x == 1){
		    print OUT ";";
		  }
		  else{
		  print OUT "$temp[0]_$temp[1]_$temp[3]/$temp[4]\t";
		  }
		  print OUT "$a=$out{$a}";
		  $x=1;
		}
	      }
	      if ($x == 1){print OUT "\n";}
	    }
	  }
	}
      }
    }
    close F;
  }
}
close LS;


#add in the sites where all background is reference
my $chr='none';my @seq=(1);
open FA,("$reference");
while (<FA>){
  if (/^>(\S+)/){
    my $mychr=$1;
    foreach my $a (keys(%positions1)){
      if ($positions1{$a} == 1){
	if ($a =~ /(\S+):(\S+)/){
	  my $two=$2;my $one=$1;
	  if ($chr eq $one){
	    if ($seq[$two] =~ /[ATCG]/){
	      print EU "$one","_$two","_$seq[$two]/$seq[$two]\t";
	      print OUT "$one","_$two","_$seq[$two]/$seq[$two]\t";
	      if ($chr eq 'Y'){
		print EU "0=",$euno{$chr},"\n";
		print OUT "0=",$outno{$chr},"\n";
	      }
	      elsif($chr eq 'X'){
		print EU "00=",2*($euno{$chr}-$euno{'Y'}),",0=$euno{'Y'}","\n";
		print OUT "00=",2*($outno{$chr}-$outno{'Y'}),",0=$outno{'Y'}","\n";
	      }
	      else{
		print EU "00=",2*$euno{$chr},"\n";
		print OUT "00=",2*$outno{$chr},"\n";
	      }
	    }
	  }
	}
      }
    }
    $chr=$mychr;@seq=(1);
  }
  else{
  my @temp=split //,$_;
  unless ($temp[scalar(@temp)-1] =~ /\w/){pop @temp;}
  push @seq,@temp;
}
}

close FA;


close EU;
close OUT;
