#!/usr/bin/perl -w
use strict;
use Statistics::Zed;

#This was run first and bad files excluded, then 23andme_similarity_filter.pl was run (after 23andme_similarity_matrix.pl) and then this was run again on the cleaned files
#on the second run set the $cleanfile parameter to the output from 23andme_similarity_filter.pl
#to compare european to all use oneliner at end of this file, after running this once on each
my $cleanfile; # = '/home/gough/data/rihab/filelist_clean.list';
my $usage="23andme_similarity_filter.pl <directory with 23andme files/directories> <#list of clean files filelist_clean.list>";
die "Usage: $usage\n" unless (@ARGV >= 1);

#first list 23andme files
my $directory='/home/gough/data/rihab/23andmedata/';
my $prepend;
my $suffix;
if (defined($ARGV[1])){$cleanfile=$ARGV[1];$suffix='_clean';}else{$suffix='';}
if (defined($ARGV[0])){$directory=$ARGV[0];$prepend=$directory;}else{$prepend='';}
my $background='/home/gough/data/snowflake/background_euro_23andme.dat';
my $output='/home/gough/data/rihab/'.$prepend."zscores$suffix.tab";
my $output2='/home/gough/data/rihab/'.$prepend."expectation$suffix.tab";
my $badfile='/home/gough/data/rihab/'.$prepend."badfiles$suffix.tab";
if ($prepend =~ /^\/home/){
$output=$prepend."zscores$suffix.tab";$output2=$prepend."expectation$suffix.tab";$badfile=$prepend."badfiles$suffix.tab";
}
my $cutoff=35;
my (@files,@bases);

if (defined($cleanfile)){
  open C,("$cleanfile");
  while (<C>){
  if (/$directory(\S+_Full_\d+.txt)/ or /$directory(\S+\.23andme.txt)$/ or /$directory(\S+\.tsv)$/ or /$directory(\S+\.23andme)$/){
    push @files, $1;
  }
  }
  close C;
print STDERR "loaded ",scalar(@files)," 23andme files from $cleanfile\n";
}
else{
open LS,("find $directory |");
while (<LS>){
  if (/$directory(\S+_Full_\d+.txt)/ or /$directory(\S+\.23andme.txt)$/ or /$directory(\S+\.tsv)$/ or /$directory(\S+\.23andme)$/){
    push @files, $1;
  }
  else{
print STDERR "Skipping (not a 23andme file): $_";
  }
}
close LS;
print STDERR "listed 23andme files\n";
}

#load in the background allele frequencies
my %alleles2;
my %positions2;
open F,("$background");
while (<F>){
  if (/^(\S+)_(\d+)_([ATCG]+)\/([ATCG]+)\s+(\S+)/){
    my $id="$1:$2";
    my $zero=$3;
    my $one=$4;
    unless (length($zero) > 1 or length($one) > 1){
    $positions2{$id}=1; 
      my @temp=split /;/,$5;
      foreach my $type (@temp){
	if ($type =~ /^00\=(\d+)/){
	  $alleles2{$id}{$zero}+=2*$1;
	}
	elsif ($type =~ /^11\=(\d+)/){
	  $alleles2{$id}{$one}+=2*$1;
	}
	elsif ($type =~ /^01\=(\d+)/ or $type =~ /^10\=(\d+)/){
	  $alleles2{$id}{$zero}+=$1;
	  $alleles2{$id}{$one}+=$1;
	}
	elsif ($type =~ /^0\=(\d+)/){
	  $alleles2{$id}{$zero}+=$1;
	}
	elsif ($type =~ /^1\=(\d+)/){
	  $alleles2{$id}{$one}+=$1;
	}
	else{
	  unless ($type =~ /^\.\.\=\d+/ or $type =~ /^\d\d?=\d+/){
	    print STDERR "Failed to parse this type: $type\n";
	  }
	}
      }
    }
  }
  else{
print STDERR "Failed to parse this line: $_";
  }
}
close F;
print STDERR "loaded background file\n";



#load in the 23andme allele frequencies
my %alleles1;
my %positions1;
foreach my $file (@files){
print STDERR "$directory/$file\n"; 
    open F,("$directory/$file");
while (<F>){
  if (/\S/){
    unless (/^#/){
      my @columns=split  /\t/,$_;
      my $id="$columns[1]:$columns[2]";
      if (exists($positions2{$id})){
	chomp($columns[3]);
	my @letters=split //,$columns[3];
	foreach my $letter (@letters){
	  $letter=uc($letter);
	  if ($letter =~  /[ATCG]/){
	    $alleles1{$id}{"$letter"}++;
	    $positions1{$id}=1;
	  }
	}
      }
    }
  }
}
close F;
}
print STDERR "loaded ",scalar(@files)," 23andme files\n";



#go through every position
my %badbase;
my $multiallelic=0;
my $scores=0;
my @zscores;
open OUT,(">$output");
foreach my $position (keys(%positions1)){
  @bases=();
  #check alleles are simple
  if (exists($alleles1{$position}{"A"}) or exists($alleles2{$position}{"A"})){
    push @bases, 'A';
  }
  if (exists($alleles1{$position}{"T"}) or exists($alleles2{$position}{"T"})){
    push @bases, 'T';
  }
  if (exists($alleles1{$position}{"G"}) or exists($alleles2{$position}{"G"})){
    push @bases, 'G';
  }
  if (exists($alleles1{$position}{"C"}) or exists($alleles2{$position}{"C"})){
    push @bases, 'C';
  }
  if (scalar(@bases) > 2){
    $multiallelic++;
    my $breakdown=&Multis($position);
    print "#multiallelic at $position: $breakdown\n";
      print OUT "$position\t","-","\tmultiallelic\n";
  }
  elsif (scalar(@bases) > 1){
      #setting missing hash values to zero
      unless (exists($alleles1{$position}{$bases[0]})){
	$alleles1{$position}{$bases[0]}=0;
      }
      unless (exists($alleles1{$position}{$bases[1]})){
	$alleles1{$position}{$bases[1]}=0;
      }
      unless (exists($alleles2{$position}{$bases[0]})){
	$alleles2{$position}{$bases[0]}=0;
      }
      unless (exists($alleles2{$position}{$bases[1]})){
	$alleles2{$position}{$bases[1]}=0;
      }
      
      #calculating the Z statistic
      my $n1=$alleles1{$position}{$bases[0]}+$alleles1{$position}{$bases[1]};
if ($n1 == 0){next;}
      my $p1=$alleles1{$position}{$bases[0]}/$n1;
      my $n2=$alleles2{$position}{$bases[0]}+$alleles2{$position}{$bases[1]};
if ($n2 == 0){next;}
      my $p2=$alleles2{$position}{$bases[0]}/$n2;
      my $p=($alleles1{$position}{$bases[0]}+$alleles2{$position}{$bases[0]})/($n1+$n2);

      #Z score calculation
      my $Z=($p1-$p2)/((($p*(1-$p))*(1/$n1+1/$n2))**0.5);
      
      print OUT "$position\t",abs($Z),"\t$bases[0]:$alleles1{$position}{$bases[0]}/$alleles2{$position}{$bases[0]},$bases[1]:$alleles1{$position}{$bases[1]}/$alleles2{$position}{$bases[1]}\n";
      push @zscores,abs($Z);
      $scores++;

      #Store which allele is offending
      if (abs($Z) >= $cutoff){
	if ($alleles1{$position}{$bases[1]} == 0){
	  $badbase{$position}=$bases[0];
	}
	elsif ($alleles2{$position}{$bases[1]} == 0){
	  $badbase{$position}=$bases[1];
	}
	elsif ($alleles1{$position}{$bases[0]}/$alleles1{$position}{$bases[1]} > $alleles2{$position}{$bases[0]}/$alleles2{$position}{$bases[1]}){
	  $badbase{$position}=$bases[0];
	}
	else{
	  $badbase{$position}=$bases[1];
	}
#      print OUT "Check: previous line's badbase is $badbase{$position}\n\n";
      }

  }
}
close OUT;
print STDERR "There were $multiallelic multi-allelic cases found and $scores Z-scores calculated\n";


#plot the expectation versus observed for Z-values
my $zed = Statistics::Zed->new();
my $bins=1000;
my %pvalues;
my $cumulative=0;
open OUT2,(">$output2");
#$p_value = $zed->p_value($z);#check for zero Pvalue

foreach my $Z (@zscores){
  if ($Z > 37){
    $Z=37;#38 is the limit for P-value calculation in the module (returns zero)
  }
  my $bin=int($bins*$Z/37);
  $pvalues{$bin}++;
}
for my $bin (0 .. $bins){
  $bin=$bins-$bin;
  my $Z=37*$bin/$bins;
  if (exists($pvalues{$bin})){
    $cumulative+=$pvalues{$bin};
  }
  my $pvalue = $zed->p_value(value => $Z, tails => 1);

# print OUT2 $Z,"\t",2*$pvalue*$scores/$cumulative,"\n";#ratio of expected to observed vs Z score
  print OUT2 2*$pvalue*$scores,"\t$cumulative\t$Z\n";#expected vs observed Z score
}

close OUT2;

#See how many bad alleles in each file at teh set cutoff
open OUT3,(">te");
foreach my $file (@files){
  my $bads=0;
  open F,("$directory/$file");
  while (<F>){
    unless (/^#/){
      my @columns=split  /\t/,$_;
      my $id="$columns[1]:$columns[2]";
      if (exists($badbase{$id})){
	chomp($columns[3]);
	my @letters=split //,$columns[3];
	foreach my $letter (@letters){
	  $letter=uc($letter);
	  if ($letter =~  /[ATCG]/){
	    if ($letter eq $badbase{$id}){
	      $bads++;
	      last;
	    }
	  }
	}
      } 
    }
  }
  close F;
  print OUT3 "$bads\t$file\n";
}
close OUT3;
system ("sort -nr te > $badfile; rm te");
print STDERR "checked all 23andme files for bad alleles\n";



sub Multis{
  my $alls='';
  my $position = shift;
  
foreach my $base (@bases){
  if (exists($alleles2{$position}{$base})){
$alls=$alls."background $base: $alleles2{$position}{$base} ";
  }
  if (exists($alleles1{$position}{$base})){
$alls=$alls."23andme $base: $alleles1{$position}{$base} ";
  }
}

return ($alls);
}


#perl -ne '$cutoff=3;if (/^(\S+)\t(\S+)\t(\S+)/){$one=$1;$two=$2;$three=$3;if ($two eq "-"){$s{$one}=9;$t{$one}="$two\t$three";}elsif ($two >= $cutoff){if ($three =~ /\w\:(\d+)\/(\d+)\,\w\:(\d+)\/(\d+)/){if ($1/($3+0.1) > $2/($4+0.1)){$s{$one}=1;}else{$s{$one}=0;}$t{$one}="$two\t$three";}else{print $_;}}}END{open F,("zscores_clean.tab");while (<F>){if (/^(\S+)\t(\S+)\t(\S+)/){$one=$1;$two=$2;$three=$3;if ($two eq "-"){$ss=9;}elsif ($two >= $cutoff){if ($three =~ /\w\:(\d+)\/(\d+)\,\w\:(\d+)\/(\d+)/){if ($1/($3+0.1) > $2/($4+0.1)){$ss=1;}else{$ss=0;}if ((exists($s{$one}) and $ss == $s{$one}) or $ss == 9 or $s{$one} == 9){print "$one\t$two\t$three\t$t{$one}\n";}}else{print STDERR "not parsed:\t",$_;}}}}close F;}' all/zscores_clean.tab > bad_alleles.tab 
