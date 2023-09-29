#!/usr/bin/perl -w
use strict;

my (@files,@summaries);
my (%store,%lines,%pvals,%pnums,%r,%t);
my $top40=1;
my $combine=0;
my $correlates=1;
my $size=1000;
my (%prate,%trate);
my $covj;

open LS, ("ls recall$size*.* |");
while (<LS>){
if (/(\S+\.\d+)$/){push @files, $1;}
}
close LS;

#correlates and T-tests
if ($correlates == 1){
my (%cohort,%scores);my @phenos;

open F,("answers34.tsv");
while (<F>){
  my @temp=split /\t/,$_;
  if ($temp[0] =~ /^\d+$/){
    push @phenos,$temp[1];
    if (exists($cohort{$temp[1]})){$cohort{$temp[1]}.=",$temp[2]";}else{$cohort{$temp[1]}=$temp[2];}
    $scores{$temp[1]}{$temp[2]}=$temp[5];
  }
}
close F;

  my (%covars);
  open F,("covariates_raw.tsv");
  while (<F>){
    unless (/^#/){
      $covj=1;
      unless (/^\d+\t/){die "parse error for covariates_raw.tsv : $_\n";}
    my @temp = split /\t/,$_;chomp($temp[scalar(@temp)-1]);
    if ($temp[1] == 1){
      $covars{$temp[3]}{$covj}=1;#sex 1
      $covars{$temp[3]}{$covj+1}=0;#sex 2
    }
    elsif($temp[1] == 2){
      $covars{$temp[3]}{$covj}=0;#sex 1
      $covars{$temp[3]}{$covj+1}=1;#sex 2
    }
    else{
      $covars{$temp[3]}{$covj}=0;
      $covars{$temp[3]}{$covj+1}=0;
    }
    $covj++;
foreach my $i (4 .. 13){
  $covj++;
  $covars{$temp[3]}{$covj}=$temp[$i];
  $covj++;
  $covars{$temp[3]}{$covj}=-$temp[$i];
}
      $covj++;
my $x=int($temp[2]*7+0.0001);
    for my $i (0 .. 5){
      if (int($temp[2]*7+0.0001) == $i){$covars{$temp[3]}{$covj+$i}=1;}#array type 1-6
      else{$covars{$temp[3]}{$covj+$i}=0;}
    }
      $covj+=5;
  }
  }
close F;$covj-=7;#removing the array type covariates - make it -6 if you can fix the missing last PC10 that was overwritten (covj = 22) by bug

foreach my $pheno (@phenos){$trate{$pheno}=0;
  my @peeps=split /,/,$cohort{$pheno};
  for my $cov (1 .. $covj){
    my ($n,$ym,$xm,$rt,$rb1,$rb2,$x,$y,$t,$r);
    foreach my $id (@peeps){
    $n++;$xm+=$scores{$pheno}{$id};$ym+=$covars{$id}{$cov};
  }
    $ym=$ym/$n;$xm=$xm/$n;
   foreach my $id (@peeps){
$x=$scores{$pheno}{$id};$y=$covars{$id}{$cov};$rt+=($x*$y);$rb1+=$x**2;$rb2=$y**2;
}
    if (($rt-$n*$xm*$ym) == 0){
$r{$pheno}{$cov}=0;
    }
    else{
      $r{$pheno}{$cov}=($rt-$n*$xm*$ym)/(($rb1-$n*$xm**2)-($rb2-$n*$ym**2))**0.5;
    }
    $t{$pheno}{$cov}=$r{$pheno}{$cov}*(($n-2)/(1-$r{$pheno}{$cov}**2))**0.5;
     if ($r{$pheno}{$cov} >= 0.05){$prate{$cov}+=1/scalar(@phenos);$trate{$pheno}=1;}
}
}

}



if ($combine == 1 or $top40 == 1){
  #load data
my $first=1;
foreach my $file (@files){
  open F,("$file");
  while (<F>){
    unless (/^#/){
      my @temp=split /\t/,$_;
      if (scalar(@temp >= 9)){
      if  ($temp[0] =~ /(\S+)/){$temp[0]=$1;}else{print STDERR "line bad: @temp\n",$_;}
      my $id = "$temp[0]|$temp[2]";
      if ($first == 1){unless (exists($lines{$temp[0]})){push @summaries,$temp[0];}$lines{$temp[0]}++;}
      $store{$id}{1}+=$temp[1]/scalar(@files);
      $store{$id}{3}+=$temp[3]/scalar(@files);
      $store{$id}{5}+=$temp[5]/scalar(@files);
      $store{$id}{8}+=$temp[8]/scalar(@files);
      $pvals{$temp[7]}{$temp[0]}+=$temp[1];$pnums{$temp[7]}{$temp[0]}++;
 if (exists($store{$id}{6})){$store{$id}{6}.=",$temp[5]|$temp[6]";}else{$store{$id}{6}="$temp[5]|$temp[6]";}
      if (exists($store{$id}{7})){
	my @t=split /,/,$store{$id}{7};
	my $done=0;
	for my $i (0 .. scalar(@t)-1){
	  my ($one,$two);
	 if ( $t[$i]=~/^(\S+)\|(\d+)/){$one=$1;$two=$2;}else{print STDERR "problem with t: -$t[$i]-   store: -$store{$id}{7}-\n";}
	  if ($temp[7] eq $one){
	    $done=1;$two++;
	    $t[$i]="$one|$two";
	  }
	}
	if ($done == 1){$store{$id}{7} = join ',',@t;}else{$store{$id}{7}="$store{$id}{7},$temp[7]|1";}
      }
      else{
	$store{$id}{7}="$temp[7]|1";
      }
    }
    }
  }
  close F;
  $first=0;
}

#make averages
foreach my $s (@summaries){
  my $maxZ=0;my $Z;
  for my $i (1 .. $lines{$s}){
    my $id = "$s|$i";unless (exists($store{$id}{6})){print STDERR "missing store 6: -$id-\n";}
    my $var=0;
    my @temp=split /,/,$store{$id}{6};
    foreach my $l (@temp){
      if (  $l =~ /(\S+)\|(\S+)/){
	$var+= (($2)**2+($1-$store{$id}{5})**2)/(scalar(@temp));
      }
      else{
print STDERR "l is wrong:  -$l-   id: -$id-  store:  $store{$id}{6}\n";
      }
    }
    $store{$id}{6}=$var**0.5;
    if ($store{$id}{6} == 0){$store{$id}{4}=0;}
    else{$store{$id}{4}=($i-$store{$id}{5})/$store{$id}{6};}
    if (($store{$id}{4})**2 > $maxZ**2){$maxZ=$store{$id}{4};$Z=$id;}
  }
    print STDERR "$Z\t$maxZ\t$store{$Z}{5}\n";
}


#output data
if ($combine == 1){
  foreach my $s (@summaries){
  for my $i (1 .. $lines{$s}){
print "$s\t";
my $id = "$s|$i";
print "$store{$id}{1}\t$i";
for my $j (3 .. 6){
print  "\t$store{$id}{$j}";
}
print "\t";
my @temp=split /,/,$store{$id}{7};
my $first=1;my $most=0;my $mostcount=0;
foreach my $i (@temp){
  if ($i =~ /(\S+)\|(\S+)/){
    if ($2 > $mostcount){$mostcount=$2;$most=$i;}
    if ($2 > scalar(@files)/10){
      unless ($first == 1){print ",";}
      print "$i";
      $first=0;
    }
  }
}
if ($first == 1){print $most;}
print "\t$store{$id}{8}\n";
  }
}
}
}

if ($top40 == 1){
#get top 40 nomaly phenotypes from DTC
my (@phenos);
open F,("dtc_goodones.tsv");
while (<F>){
  my @temp = split /\t/,$_;
  unless ($temp[7] eq 'term'){
    push @phenos,$temp[7];
  }
}
close F;
my $c=0;
foreach my $p (@phenos){
  $c++;if ($c == 41){print "\n\n";}
  print  "$p";my $j=0;
  foreach my $sum (@summaries){
    $j++;
    unless ($j > $covj+1){
    my $val='-';    my $z='-';
    if (exists($pnums{$p}{$sum})){
      $val=$pvals{$p}{$sum}/$pnums{$p}{$sum};
    
   

    #get Zscore
    open F,("recall$size.all");
    while (<F>){
      my @temp = split /\t/,$_;
      if ($sum eq $temp[0]){
      if ($temp[1] <= $val){
$z=$temp[4];
      }
    }
    }
      close F;
    }
    
    
      if ($sum =~ /Summary3:(\d+)/){
    if ($val <= 0.05 and $z >= 1.65){
#    if ($r{$p}{$1} >= 0.05 or ($val <= 0.05 and $z >= 1)){#if you change this also change upstairs
      print  "\t$1:$val:$z:$r{$p}{$1}:$t{$p}{$1}";
    }
        else{
print  "\t      ";
    }
      }
      else{
	print  "\t$val";
      }

  }
}
  print  "\n";
}


print STDERR "\nOverall phenotype rate: ";my $r=0;
foreach my $p (keys(%trate)){if ($trate{$p} == 1){$r++;}}
print STDERR $r/scalar(keys(%trate)),"\n";
print STDERR "Covariate rates:";
foreach my $cov (@summaries){
      if ($cov =~ /Summary3:(\d+)/){
	$cov=$1;
	unless (exists($prate{$cov})){$prate{$cov}=0;}
	print STDERR "\t$prate{$cov}";
      }
    }
print STDERR "\n";



}

