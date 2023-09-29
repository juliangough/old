#!/usr/bin/perl -w
use strict;

my @files;
my (%everything,%alleles,%badfiles,%totalcount,%totalcompared);
my $x=0;
my ($badsnps,$snps);
my $tempdir='similarity';unless (-e $tempdir){system("mkdir $tempdir");}
my $cpus=70;
my $outfile='similarity_matrix.tab';
my $outfile2='max_similarity_matrix.tab';
my $skiptempfiles=0;
my $ancestrydna=1;  #set to 1 for ancestryDNA or 0 for 23andme


open L,("filelist.list");
while (<L>){
  if (/^\d+\t(\S.*)$/){
    my $file=$1;chomp($file);push @files,$file;
}
else{
print STDERR "Failed to parse file: $_";
}
}
close L;

open S,("/data/snow/SNV/daemon/allcluster.tab");
while (<S>){
if (/^(\S+)/){
    $alleles{$1}=1;
}
}
close S;


foreach my $file (@files){
  if (-e "$file"){
 $snps=0;
open F,("$file");
while (<F>){
 unless (/^#/){
   my @temp=split /\t/,$_;chomp($temp[3]);
  if (exists($alleles{"$temp[1]:$temp[2]"})){
   $everything{$file}{"$temp[1]:$temp[2]"}=$temp[3];
}
 }
}
close F;
$x++;
  print STDERR "$x of ",scalar(@files)," files read in snps $file\n";
}
}

print STDERR "\n",scalar(keys(%alleles))," alleles mapped\n";

my $i=0;
my $blocksize=int(scalar(keys(%alleles))/$cpus+1);
foreach my $allele (keys(%alleles)){
  if ($i+1 > scalar(keys(%alleles))){
close F;
    last;
  }
  if ($i/$blocksize == int($i/$blocksize)){
    unless ($i == 0){
      close F;
    }
      my $tempfile="$tempdir/kfhdhfhkjmlociwj.".$i/$blocksize;
open F,(">$tempfile");
  }
  $i++;
  foreach my $file (@files){
    if (exists($everything{$file}{$allele})){
      print F "$everything{$file}{$allele},";
    }
    else{
print F ",";
    }
  }
  print F "$allele\n";
}

%everything=();
%alleles=();

print STDERR "Temporary files created\n";


#-------FORK-----------------------------------------------
my ($pid,$childid);
foreach my $cpu (0..$cpus-1){
  if (!defined($pid=fork())){
die "Fork Error!";
  }
elsif ($pid==0){
&Output($cpu);
exit 0;
  }
else{
$childid=$pid;
}
#child never gets here
}

while(wait() > 0){
;
}
#-----------------------------------------------------------


for my $output (0 .. $cpus-1){
  open F,("$tempdir/kfhdhfhkjmlociwj.$output.out");
  while (<F>){
    if (/^(\S+)\t(\S+)\t(\S+)\t(\S+)/){
$totalcount{$1}{$2}+=$3;
$totalcompared{$1}{$2}+=$4;
    }
  }
  close F;
}

open OUT,(">$outfile");
open OU2T,(">$outfile2");
foreach my $i (0 .. scalar(@files)-2){
  my $max=0;
     foreach my $j ($i+1 .. scalar(@files)-1){
      unless (exists($totalcompared{$i}{$j})){
	$totalcompared{$i}{$j}=0;
      }
       if ($totalcompared{$i}{$j} == 0){
	 print OUT "$i\t$j\t0\n";
       }
       else{
	 print OUT "$i\t$j\t",$totalcount{$i}{$j}/$totalcompared{$i}{$j},"\n";
	 if ($totalcount{$i}{$j}/$totalcompared{$i}{$j} > $max){
	   $max=$totalcount{$i}{$j}/$totalcompared{$i}{$j};
	 }
       }
    }
  print OU2T "$max\t$i\n";
 }
close OU2T;
close OUT;

system ("sort -n $outfile2 > te ; mv te $outfile2");
system ("rm -rf $tempdir");

sub Output{
  my $id = shift;
  my $file="$tempdir/kfhdhfhkjmlociwj.$id";
  my %count;
  my %compared;

  unless (-e "$file.out"){
  open F,("$file");
  while (<F>){
    my @temp=split /,/,$_;
    pop(@temp);
    foreach my $i (0 .. scalar(@temp)-2){
      foreach my $j ($i+1 .. scalar(@temp)-1){
	unless ($temp[$i] eq '' or $temp[$j] eq '' or $temp[$i] eq '--' or $temp[$j] eq '--'){
	  $compared{$i}{$j}++;
	  my @tempi=split //,$temp[$i];my @tempj=split //,$temp[$j];
	  my $common=0;
	  foreach my $ti (@tempi){
	    foreach my $tj (@tempj){
	      if ($tj eq $ti){
		$common++;
	      }
	    }
	  }

	  if ($common == 4){#homozygous match
	    $count{$i}{$j}++;
	  }
	  elsif($common == 1 and scalar(@tempi) == 1 and scalar(@tempj) == 1){#single allele match
	    $count{$i}{$j}++;
	  }
	  elsif (($common == 1 or $common == 2) and (scalar(@tempi) == 1 or scalar(@tempj) == 1)){#a single and a double allele with a homozygous or heterozygous match
	    $count{$i}{$j}+=0.5;
	  }
	  elsif($common == 2 and $tempi[0] ne $tempi[1] and $tempj[0] ne $tempj[1]){#heterozygous match
	    $count{$i}{$j}++;
	  }
	  elsif ($common == 2){#heterozygous and homozygous
	    $count{$i}{$j}+=0.5;
	  }
	}
      }
    }
  }
  close F;
  open G,(">$file.out");
  foreach my $i (0 .. scalar(@files)-2){
    foreach my $j ($i+1 .. scalar(@files)-1){
      unless (exists($count{$i}{$j})){
	$count{$i}{$j}=0;
      }
      unless (exists($compared{$i}{$j})){
	$compared{$i}{$j}=0;
      }
       print G "$i\t$j\t$count{$i}{$j}\t$compared{$i}{$j}\n";
    }
  }
  close G;
}
}


