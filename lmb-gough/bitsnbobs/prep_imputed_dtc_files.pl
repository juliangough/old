#!/usr/bin/perl -w
use strict;

my $num;my $chromosome;
if (defined($ARGV[0])){$num=$ARGV[0];}else{$num=-1;}
#JulianGough Sat 15th Oct 2022 in Cafe Nero Cambridge, Regent Street whilst Juno is at ballet

#luca files
my $btches='/home/clu/unrelated_batch.tsv';
my %leftover;

#read in the sexes
my %sex;
open S,("/data2b/snow/impute/max/impute/final/X_all.gen.samples");
while (<S>){
  if (/^(\d\d\d\d+)\s+\d+\s+\d\s+(\d)$/){
    $sex{$1}=$2;
  }
}
close S;

#read in the imputed genotypes
my $head;my %predata;my %allele;my %include;
open LS,("ls /data2b/snow/impute/max/impute/imputed_*_final.vcf |");
while (<LS>){
  my $file = $_;chomp($file);my @ids;
  open F,("$file");
  while (<F>){
    if (/^(#CHROM\s+.+\sFORMAT)\s+(\S.+\S)$/){
      $head="##fileformat=VCFv4.2\n$1";
      @ids=split /\s+/,$2;
    }
    elsif(/^(\S+)\t(\d+)(\t\S+\t)(\S+)\t(\S+)(\t\S+\t\S+\t\S+\t\S+)\t(\S.+\S)$/){
      my $loc="$1:$2:$4:$5";my $chr=$1;
      $predata{$loc}="$1\t$2$3$4\t$5$6";$include{$loc}=1;
      my @t=split /\t/,$7;
      foreach my $i (0 .. scalar(@ids)-1){
	$leftover{$ids[$i]}=1;
	if ($chr eq 'X'){
	  if ($sex{$ids[$i]} == 1){
	    if ($t[$i] =~ /^(\d)\|(\d)$/){
	      if ($1 == $2){
		$allele{$ids[$i]}{$loc}=$1;
	      }
	      else{
		print STDERR "8\t$t[$i]\n";die;
	      }
	    }
	    else{
	      if ($t[$i] eq "."){
		$allele{$ids[$i]}{$loc}=$t[$i];
	      }
	      else{
		print STDERR "9\t$t[$i]\n";die;
	      }
	    }
	  }
	  else{
	    $allele{$ids[$i]}{$loc}=$t[$i];
	  }
	}
	else{
	  $allele{$ids[$i]}{$loc}=$t[$i];
	}
      }
    }
    else{
      unless (/^#/){print STDERR "PARSE FAILED VCF: $_";}
    }
  }
  close F;
}
close LS;

#check and add missing original genotypes including Y chromosome and mitochondria
open LS,("ls /data2b/snow/impute/final/checked/*.vcf |");
while (<LS>){
  my $file = $_;chomp($file);my $id;
  my %once;
  open F,("$file");
  while (<F>){
    if (/^(#CHROM\s+.+\sFORMAT)\s+(\d+)$/){
      $id=$2;
    }
    elsif(/^(\S+)\t(\d+)(\t\S+\t)(\S+)\t(\S+)(\t\S+\t\S+\t\S+\t\S+)\t(\S+)$/){
      my $loc="$1:$2:$4:$5";my $chr=$1;
      unless (exists($predata{$loc})){$predata{$loc}="$1\t$2$3$4\t$5$6";}
      my $t=$7;

      if (exists($allele{$id}{$loc})){
	unless ($allele{$id}{$loc} eq $t or $t eq '.' or ($t eq '1|0' and $allele{$id}{$loc} eq '0|1') or ($t eq '0|1' and $allele{$id}{$loc} eq '1|0')){
	  if ($allele{$id}{$loc} eq '.'){
	    $allele{$id}{$loc}=$t;
	  }
	  else{
	    unless ($t eq "$allele{$id}{$loc}|$allele{$id}{$loc}"){
print STDERR "disagree: orig: $t imputed: $allele{$id}{$loc} loc: $loc id: $id\n";
$allele{$id}{$loc}=$t;
}
	  }
	}
      }
      else{
	unless ($t eq '.'){
	$allele{$id}{$loc}=$t;
	$include{$loc}=1;
      }
      }
    }
    else{
      unless (/^#/){print STDERR "PARSE FAILED VCF: $_";}
    }
  }
  close F;
}
close LS;



for my $batch (0 .. 5){
  if ($num > -1){$batch = $num;}
#grab the list for this batch
my @list;
open B,("$btches");
while (<B>){
  if (/^(\d+)\s+(\d)$/){
    if ($2 == $batch){
      push @list,$1;
    }
  }
  else{
print STDERR "UNPARSED: $_";
  }
}
close B;
#----------------------------

  if ($batch == 5){
    unless (scalar(@list) == 0){die;}
open B,("$btches");
while (<B>){
  if (/^(\d+)\s+(\d)$/){
    $leftover{$1}=0;
  }
}
close B;
foreach my $i (keys(%leftover)){
  if ($leftover{$i} == 1){
      push @list,$i;
  }
}
  }

unless (-e "/data2b/snow/impute/final/batch$batch"){  system ("mkdir /data2b/snow/impute/final/batch$batch");}
  foreach my $id (@list){
  open OUT, (">/data2b/snow/impute/final/batch$batch/$id.vcf");
  print OUT $head;
print OUT "\t$id";
  print OUT "\n";
  foreach my $loc (keys(%include)){
    my $line = $predata{$loc};my $blank=0;
      if (exists($allele{$id}{$loc})){
	$line .= "\t$allele{$id}{$loc}";
	unless ($allele{$id}{$loc} eq '.'){$blank=1;}
      }
      else{
$line .= "\t.";
	}
    $line .= "\n";
    unless  ($blank == 0){
      print OUT $line;
    }
  }
  close OUT;
}
  if ($num > -1){last;}
}
