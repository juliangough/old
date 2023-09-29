#!/usr/bin/perl -w
use strict;

my $num=1;my $chromosome;
if (defined($ARGV[0])){$num=$ARGV[0];}else{$num=0;}
if (defined($ARGV[1])){$chromosome=$ARGV[1];}else{$chromosome=0;}
#JulianGough Sat 8th Oct 2022 in Cafe Nero Cambridge, Regent Street whilst Juno is at ballet
my $exonfilter=1;#set to filter exons or not - takes too long
print STDERR "WARNING: does not print lines where every entry is the same, including where all at heterozygous or homozygous rare!\n\n";
#my %s;
#open L,("ls /data/gough/imputation/*_1.phased.impute2 |");
#while (<L>){
#my $file =$_;chomp($file);
#open F,("$file");
#while(<F>){
#if (/^\S+\s+\S+\s+(\d+)\s+(\S+)\s+(\S+)\s+(\S.*\S)$/){
#$s{$chr}{$1}=0;
#}
#}
#close F;
#}
#close L;


# #read exon regions
# my %exon;my %next;
# if ($exonfilter == 1){
# my $chr='';my $prev=0;
# open EXONS,("/home/snow/data/imputation/hg19_exons.all.bed");
# while (<EXONS>){
#   if (/^NC_0000(\d\d)\.\S+\s+(\d+)\s+(\d+)\s/){
# my $c=$1;my $a=$2;my $b=$3;
# if ($c =~ /0(\d)/){$c=$1;}
#  if ($c == 23){$c='X';}elsif($c == 24){$c='Y';}
# if ($chr ne $c){unless ($chr eq ''){print STDERR "Chr: $chr\n";}$prev=0;}
#     $chr=$c;
# $exon{$chr}{$a}=$b;
# unless ($prev == 0){$next{$chr}{$prev}=$a;}$prev=$b;
#   }
#   else{
# #print STDERR $_;
# }
# }
# close EXONS;
# }

my (%first,%span,%next);
my $last;
open EXONS,("/beegfs3/gough/hg19_exons.all.bed");
while (<EXONS>){if (/^NC_0+(\d+)\.\d+\s+(\d+)\s+(\d+)/){
  unless (exists($first{$1})){
    $first{$1}=$2;$span{$1}{$2}=$3;
    $last=$3;
  }
  elsif(exists($span{$1}{$2})){
    if ($span{$1}{$2} == $3){
      print STDERR "WARNING: $2 $3 $last repeat line\n";
      next;
    }
    else{
      print STDERR "out of sequence and not repeated line: $2 $3 $last $_";
      if ($span{$1}{$2} < $3){
	$span{$1}{$2}=$3;
	$last=$3;
      }
      next;
    }
  }
  elsif ($2 < $last){
    print STDERR "WARNING: out of seqeunce line: $last $1 $2 $3 $_";
  }
  else{
    $span{$1}{$2}=$3;
    $next{$1}{$last}=$2;
    $last=$3;
  }
}
	  }
close EXONS;

#add allele codes
my %allele = ('100','0|0','010','1|0','001','1|1');

#do the thing
for my $chunk  (1 .. 9){
if ($num > 0 and $num < 10){$chunk=$num;}
my @files;
if ($chromosome ne '0'){
  open VCF, (">/beegfs3/gough/imputation/dtc/impute/imputed_$chromosome\_final.vcf");
}
else{
  open VCF, (">/home/snow/data/imputation/dtc/impute/imputed_final_$chunk.vcf");
}
#get headers for VCF files
open F,("/beegfs3/gough/imputation/dtc/impute/final/$chromosome\_all.vcf");
while (<F>){
if (/^\#/){
print VCF $_;
}
}
close F;

#get imputed content
my %ends;my %map;
if ($chromosome ne '0'){
  open LS,("ls /beegfs3/gough/imputation/dtc/impute/imputed/$chromosome\_*.phased.impute2 |");
}
else{
  open LS,("ls /home/snow/data/imputation/dtc/impute/imputed/*_$chunk.phased.impute2 |");
}
while (<LS>){
if (/\/imputation\/dtc\/impute\/imputed\/(\S+)_(\d+)_(\d+).phased.impute2$/){
my $w=$1;my $x=$2;my $z=$3;
if ($w eq 'X'){$w = 23;}if ($w eq 'Y'){$w = 24;}
until (length($x) == 10){$x='0'.$x;}
#until (length($w) == 10){$w='0'.$w;}
$map{"$w$x"}=$_;chomp($map{"$w$x"});
$ends{"$w$x"}=$3;
}
}
close LS;

my @t=keys(%map);
my @temp = sort  numerically @t;
foreach my $f (@temp){
push @files,$map{$f};
}
print STDERR scalar(@files)," files for chunk $chunk\n";

foreach my $file (@files){

if ($file =~ /\/imputation\/dtc\/impute\/imputed\/(\S+)_(\d+)_(\d+).phased.impute2$/){
  my $chrome=$1;my $start=$2;my $stop=$3;my $chr=$chrome;
  if ($chr eq "X"){$chr=23;}elsif($chr eq "Y"){$chr=24;}
#my $end=$exon{$chr}{$start};
open GEN,("$file");
while (<GEN>){
if (/^\S+\s+(\S+)\s+(\d+)\s+(\S+)\s+(\S+)\s+(\S.*\S)$/){

my $rsid=$1;
my $pos=$2;
my $maj=$3;
my $min=$4;
my @temp = split /\s+/,$5;
if ($rsid =~ /^(\S+):\d+:\S+:\S+$/){$rsid=$1;}
#print STDERR "$chr\t$start\t$stop\t$pos\t$end\t$next{$chr}{$end}\t$exon{$chr}{$next{$chr}{$end}}\n";
#my $go=0;if($pos >= $next{$chr}{$end}){until ($pos <= $end){$end=$exon{$chr}{$next{$chr}{$end}};}}if ($pos <= $end){$go=1;}
#my $go=1;

my $out=0;my $go=0;
if ($pos < $first{$chr}){
  $out=1;
}
my $a=$first{$chr};
until ($out == 1){
#  if ($pos == 19164633){print STDERR "debug: $chr $a $out $flag $span{$chr}{$a} $next{$chr}{$span{$chr}{$a}} $first{$chr} $pos\n";}
  if ($pos < $span{$chr}{$a}){
    $go=1;
    $out=1;
  }
  elsif (exists($next{$chr}{$span{$chr}{$a}})){
    if ($pos > $next{$chr}{$span{$chr}{$a}}){
      $a=$next{$chr}{$span{$chr}{$a}};
    }
    else{$out=1;
       }
  }
  else{
    $out=1;
  }
}


if ($go == 1){
my $line = "$chrome\t$pos\t";
if (length($rsid) > 3 ){$line .=  "$rsid";}else{$line .= ".";}
$line .= "\t$min\t$maj\t.\t.\t.\tGT";
unless (scalar(@temp)/3 == int(scalar(@temp)/3)){print STDERR "FAIL ALLELE PARSE: ","@temp,","\n";}
my $nonzero=0;my $last='s';
for my $j (0 .. (scalar(@temp)/3)-1){
my $code='';
for my $k (0 .. 2){
if ($temp[($j*3+$k)] >= 0.994){
$code.='1';
}
elsif($temp[($j*3+$k)] <= 0.006){
$code.='0';
}
}
if (exists($allele{$code})){
if ($allele{$code} ne $last and $last ne 's' ){$nonzero = 1;}$last=$allele{$code};
$line.= "\t$allele{$code}";
}
else{
$line.= "\t.";
}
}
if ($nonzero == 1){print VCF "$line\n";}
}
}
else{
print STDERR "FAIL GEN PARSE $file:$_";
}
}
close GEN;
}
}
close VCF;
if ($num > 0 and $num < 10){last;}
}


#SORT-BY-NUMBER---------------------------------------
sub numerically {$a <=> $b;}
#-----------------------------------------------------
