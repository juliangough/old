#!/usr/bin/perl -w
use strict;
#use Statistics;
use List::Util qw(shuffle);


my ($arg1,$arg2);

my $data='/Users/gough/DropBox/Projects/clu/questions_hasyes_recall_bootstrap.tsv';$data='';#don't use this file again
my $answers='/Users/gough/DropBox/Projects/clu/answers34.tsv';
#$answers='/data2/gough/clu/answers34.tsv';
#$answers='/data2/gough/clu/ddd_answers_close.tsv';
#$answers='/data2/gough/clu/ddd_answers_exact.tsv';
#$answers='/data2/gough/clu/julian_ddd/answers_1133_exact.tsv';
#$answers='/data2/gough/clu/julian_ddd/answers_4239_exact.tsv';
if (defined($ARGV[0])){$answers=$ARGV[0];}
my %bin;
my @temp;
my $bins=20;
my $scorebinsize=10;
my $asymmetry=0;
my $aveperm=0;
my $bootstrapplot='no';
my $scatterplot = 'no';
my $thresholdplot = 'no';
my $decoystats='no';
my $answersscatterplot='no';
my $truerate='no';
my $scoretruerate='no';
my $decoyscatter = 'no';
my $recallyesrates  = 'no';
my $permutation_test_global  = 'no';
my $permutation_test_question  = 'no';
my $permutation_test_person  = 'no';
my $yesrate_byscorethreshold = 'no';
my $permutation_bytopquestions = 'yes';
my $permutation_byquestion_top_5 = 'no';
my $permutation_test_decoyszero = 'no';
my $top=((100)/(130.25))*0.05;
my $perms = 1000;
my $zero=0.01;
my $threshold = 0.022;
my $ratefilter=0.05;
my %termyesrates=();
my %termyesrates_norecall=();
my $type='a';
if ($type eq 'f' or $type eq 'e'){&GetDecoyYesRates();}
else{&GetYesRates();}
if (defined($ARGV[0])){$arg1=$ARGV[0];}else{$arg1=1;}
if (defined($ARGV[1])){$arg2=$ARGV[1];}else{$arg2=50}



if ($permutation_bytopquestions eq 'yes'){
%termyesrates=();&GetDecoyYesRates();
  my (%term,%isyes,%score,%decoy,%sums);
open F,("$answers");
while (<F>){
  unless (/^\S*\s+ont/){
    @temp=split /\t/,$_;chomp($temp[scalar(@temp)-1]);
          unless (exists($termyesrates{$temp[1]})){$termyesrates{$temp[1]}=0;}
             if ($temp[8] eq 'True'){$decoy{$temp[1]}{$temp[2]}=0;}else{$decoy{$temp[1]}{$temp[2]}=1;}
      if (exists($term{$temp[1]})){$term{$temp[1]}=$term{$temp[1]}.",$temp[2]";}else{$term{$temp[1]}=$temp[2];}
    $score{$temp[1]}{$temp[2]}=$temp[5];
    if ($temp[4] eq 'True'){
      $isyes{$temp[1]}{$temp[2]}=1;
      if ($termyesrates{$temp[1]} <= $ratefilter){
	$sums{$temp[1]}+=$temp[5];
      }
    }
    elsif($temp[4] eq 'False'){
      $isyes{$temp[1]}{$temp[2]}=0;
    }
  }
}
   close F;
   
print STDERR "Questions: ",scalar(keys(%term)),"  Below ",$ratefilter*100,"\% yes rate with at least one 'yes': ",scalar(keys(%sums)),"\n";

#average score of yesses in permutations
  my %permsave;
  for my $p (1 .. $perms){
    foreach my $pheno (keys(%term)){
      my @peeps=split /,/,$term{$pheno};
      my @rands=shuffle(@peeps);
      my $sum=0;my $decoys=0;my $decoyyes=0;
      foreach my $a (0 .. scalar(@rands)-1){
	if ($decoy{$pheno}{$peeps[$a]} == 1){$decoys++;}
	if ($isyes{$pheno}{$rands[$a]} == 1){
	  $sum+=$score{$pheno}{$peeps[$a]};
	if ($decoy{$pheno}{$peeps[$a]} == 1){$decoyyes++;}	  
	}
      }
      if ($decoyyes <= $decoys*$ratefilter){
	$permsave{$pheno}{$p}=$sum;
      }
      }
  }


  my %pval;my %n;
    for my $p (1 .. $perms){
      foreach my $pheno (keys(%term)){
	unless (exists($permsave{$pheno}{$p})){$permsave{$pheno}{$p}=0;}unless (exists($sums{$pheno})){$sums{$pheno}=0;}
	unless ($sums{$pheno} == 0 and $permsave{$pheno}{$p} == 0){#if yes rate always > filter exclude
	  $n{$pheno}++;
	if ($permsave{$pheno}{$p} > $sums{$pheno}){
	  $pval{$pheno}++;
	}
	elsif($permsave{$pheno}{$p} == $sums{$pheno}){
	  $pval{$pheno}+=0.5;
	}
	  else{
	    unless (exists($pval{$pheno})){$pval{$pheno}=0;}#this increases 'n' when possibility of inclusion
	  }
      }
      }
    }
my $number=0;my @keys;my @nums;
       foreach my $pheno (keys(%pval)){
	 push @nums,($perms*$pval{$pheno}/$n{$pheno});#this reduces significance of pval in the case where filtered on decoy yes rate
	 push @keys,$pheno;
       }
my @sortedkeys=OrderArray(@nums);
foreach my $pheno (@sortedkeys){
  my $pvalue=$nums[$pheno];
  $pheno=$keys[$pheno];
	 $number++;
  print "Summary\t",$pvalue/$perms,"\t$number\t",scalar(keys(%pval))*$pvalue/$perms,"\t$pheno\t$pvalue\t",scalar(keys(%pval)),"\n";
}
}


#SUB-ROUTINE-----------------------------------------------
sub GetYesRates{
my (%y,%n,%ynr,%nnr);my @temp;
  #including recall
  
open F,("$answers");
while (<F>){
  unless (/^\S*\s+ont/){
    @temp=split /\t/,$_;chomp($temp[scalar(@temp)-1]);
    if ($temp[8] =~ /\S+/){
      $nnr{$temp[1]}++;
      if ($temp[4] eq 'True'){
$ynr{$temp[1]}++;
      }
    }
      $n{$temp[1]}++;
      if ($temp[4] eq 'True'){
$y{$temp[1]}++;
      }
  }
}
  close F;
foreach my $term (keys(%nnr)){
  unless ($n{$term} >= $nnr{$term}){die;}
if (exists($ynr{$term})){$termyesrates_norecall{$term}=$ynr{$term}/$nnr{$term};}else{$termyesrates_norecall{$term}=0;}
if (exists($y{$term})){$termyesrates{$term}=$y{$term}/$n{$term};}else{$termyesrates{$term}=0;}
 }


return;
}
#----------------------------------------------------------


#SUB-ROUTINE-----------------------------------------------
sub GetDecoyYesRates{
my (%y,%n,%ynr,%nnr);my @temp;
  #including recall
  
open F,("$answers");
while (<F>){
  unless (/^\S*\s+ont/){
    @temp=split /\t/,$_;chomp($temp[scalar(@temp)-1]);
    if ($temp[8] =~ /False/){
      $nnr{$temp[1]}++;
      if ($temp[4] eq 'True'){
$ynr{$temp[1]}++;
      }
    }
        unless ($temp[8] =~ /True/){
      $n{$temp[1]}++;
      if ($temp[4] eq 'True'){
$y{$temp[1]}++;
      }
    }
  }
}
  close F;
foreach my $term (keys(%nnr)){
  unless ($n{$term} >= $nnr{$term}){die;}
if (exists($ynr{$term})){$termyesrates_norecall{$term}=$ynr{$term}/$nnr{$term};}else{$termyesrates_norecall{$term}=0;}
if (exists($y{$term})){$termyesrates{$term}=$y{$term}/$n{$term};}else{$termyesrates{$term}=0;}
 }


return;
}
#----------------------------------------------------------


































if ($bootstrapplot eq 'yes'){
  #read data into bins optionally with filters
  my $n=0;
open F,("$data");
while (<F>){
  unless (/^Ontology/){
    @temp=split /\t/,$_;
    for my $i (1 .. $bins){
      if ($temp[18] <= $i/$bins and $temp[18] > ($i-1)/$bins){
	$bin{$i}++;
	$n++;
	if ($i <= $bins/2){$asymmetry++;}else{$asymmetry--;}
      }
    }
  }
}
close F;


#plot data in bins
for my $i (1 .. $bins){
  my $value=($i-0.5)/$bins;
  print "$value\t$bin{$i}\n";
}
print "\n";
for my $i (1 .. $bins){
  my $value=($i-0.5)/$bins;
  print "$value\t",$n/$bins,"\n";
}

#print out the asymmetry
print STDERR "\n$asymmetry more on the positive side out of $n\n";
}


if ($scatterplot eq 'yes'){
#plot data
open F,("$data");
while (<F>){
  unless (/^Ontology/){
    @temp=split /\t/,$_;
    my $yesrate=$temp[15]/$temp[16];
print "$yesrate\t$temp[27]\n";
  }
}
close F;
}


if ($thresholdplot eq 'yes'){
    my $n=0;
  $asymmetry = 0;
  #read data in and sort
  my (@yesrates,@perms);
open F,("$data");
while (<F>){
  unless (/^Ontology/){
    @temp=split /\t/,$_;
    push @yesrates,$temp[15]/$temp[16];
    push @perms,$temp[18];
  }
}
  close F;

  my @sorted = OrderArray(@yesrates);
  my $permcount=0;
  for my $i (@sorted){
    if ($perms[$i] <= 0.5){$asymmetry++;}else{$asymmetry--;}
    $n++;
    my $frac=$asymmetry/$n;
    $permcount+=$perms[$i];
    $aveperm=$permcount/$n;
print "$yesrates[$i]\t$aveperm\n";
}
}



if ($decoystats eq 'yes'){
  my @thresholds=(0.12,0.1,0.09,0.085,0.0825,0.08,0.075,0.07,0.065,0.0625,0.061,0.06,0.0575,0.055,0.052,0.05,0.0475,0.045,0.04,0.035,0.03,0.025,0.023,0.021,0.02,0.015,0.01,0.008,0.005,0.003,0.002,0.001,0.0001);
    for $threshold (@thresholds){
#  for my $noofqs (1 .. 25){
    print STDERR "Threshold: $threshold\n";
#    print STDERR "Number of questions: $noofqs\n";
  my ($x1,$x2,$l1,$l2,$N1,$N2);
  my ($c1,$c2);
  my %include;
  $x1=0;$x2=0;
  my %diff;
  my (%n1,%n2);
  my $id;my $lastid='';
open F,("$answers");
while (<F>){
  unless (/^\S*\s+ont/){
    @temp=split /\t/,$_;chomp($temp[scalar(@temp)-1]);
    $id=$temp[2];
    if ($lastid ne $id){
      $lastid=$temp[2];
      %include=();
      my $other=0;my $qs=0;my $dqs=0;
      open S,("/data3/snow/SNV/daemon/$id/$id.out");
      while (<S>){
	if (/^Random/){$other=1;}
	my   @temp2=split /\t/,$_;chomp($temp2[scalar(@temp2)-1]);
	if (scalar(@temp2 > 6)){
	   if ($temp2[1] >= $threshold){
#	  if ($qs < $noofqs){
   $include{$temp2[0]}=1;
   $qs++;$c1++;
 }
 elsif($other == 1 and $qs > $dqs){
   $include{$temp2[0]}=1;
$dqs++;$c2++;
 }
}
      }
      close S;
    }
    if ($temp[8] eq 'True' and exists($include{$temp[1]}) ){
      if ($temp[4] eq 'True'){
      $n1{$temp[2]}=1;
      $x1++;
      $diff{$temp[2]}++;
    }
    }
    elsif($temp[8] eq 'False' and exists($include{$temp[1]})){
      if ($temp[4] eq 'True'){
      $n1{$temp[2]}=1;
      $x2++;
      $diff{$temp[2]}--;
    }
    }
  }
}
  close F;
  #Z score
  $N1=scalar(keys(%n1));
  $N2=scalar(keys(%n1));
  #Two Poisson rate test 
  $l1=$x1/$N1;
  $l2=$x2/$N2;
  my $zscore=($l1-$l2)/($l1/$N1+$l2/$N2)**0.5;
  print STDERR "Z-score of Poisson two rate test: $zscore  $N1 $N2 $l1 $l2 $x1 $x2\n";
  # T-test
  my $sum=0;
  foreach my $id (keys(%n1)){
    $sum+=$diff{$id};
  }
    my $mean=$sum/$N1;
  my $sqrsum=0;
    foreach my $id (keys(%n1)){
    $sqrsum+=($diff{$id}-$mean)**2;
  }
  my $var=$sqrsum/$N1;
  my $t=$mean/((($var)**0.5)/($N1)**0.5);
    print STDERR "T-test of the paired samples: $t  $mean $var $c1 $c2\n";
    print "$threshold\t$t\t$zscore\t$x1\t$x2\n";
#    print "$noofqs\t$t\t$zscore\t$x1\t$x2\n";
}
}

if ($yesrate_byscorethreshold eq'yes'){
  my (@yesno,@scores,@yesrateabove,@yesrateratio);
open F,("$answers");
while (<F>){
  unless (/^\S*\s+ont/){
    @temp=split /\t/,$_;chomp($temp[scalar(@temp)-1]);
    if ($temp[8] =~ /\S+/ or 1 == 1){#includes recall
      unless (exists($termyesrates{$temp[1]})){$termyesrates{$temp[1]}=0;}
      if ($termyesrates{$temp[1]} < 0.2){
	push @scores,$temp[5];
#	if ($temp[5] > 0.1) {print STDERR $_;}
      if ($temp[4] eq 'True'){push @yesno,1;}elsif($temp[4] eq 'False'){push @yesno,0;}else{die;}
    }
  }
  }
}
  close F;

  my $n=0;my $tsum=0;
  my @sorted = OrderArray(@scores);
  foreach my $i (0 .. scalar(@sorted)-1){
    $i = $sorted[scalar(@sorted)-1-$i];
     $n++;
    if ($yesno[$i] == 1){$tsum++;}
    push @yesrateabove,($tsum/$n);
    print "$scores[$i]\t",$tsum/$n,"\n";
  }
  print "\n";
    $n=0;$tsum=0;
  foreach my $i (0 .. scalar(@sorted)-1){
    $i = $sorted[$i];
    $n++;
    if ($yesno[$i] == 1){$tsum++;}
    if ($tsum > 0){
      print "$scores[$i]\t",($yesrateabove[scalar(@sorted)-1-$n])/($tsum/$n)-1,"\n";
    }
  }
  print "\n";

  
#   @yesno=();@scores=();@yesrateabove=();@yesrateratio=();
# open F,("$answers");
# while (<F>){
#   unless (/^\S*\s+ont/){
#     @temp=split /\t/,$_;chomp($temp[scalar(@temp)-1]);
#     if ($temp[8] =~ /\S+/){#discludes recall
#       unless (exists($termyesrates_norecall{$temp[1]})){$termyesrates_norecall{$temp[1]}=0;}
#       if ($termyesrates_norecall{$temp[1]} < 0.05){
#       push @scores,$temp[5];
#       if ($temp[4] eq 'True'){push @yesno,1;}elsif($temp[4] eq 'False'){push @yesno,0;}else{die;}
#     }
#     }
#   }
# }
#   close F;

#   $n=0;$tsum=0;
#   @sorted = OrderArray(@scores);
#   foreach my $i (0 .. scalar(@sorted)-1){
#     $i = $sorted[scalar(@sorted)-1-$i];
#      $n++;
#     if ($yesno[$i] == 1){$tsum++;}
#     push @yesrateabove,($tsum/$n);
#     print "$scores[$i]\t",$tsum/$n,"\n";
#   }
#   print "\n";
#     $n=0;$tsum=0;
#   foreach my $i (0 .. scalar(@sorted)-1){
#     $i = $sorted[$i];
#     $n++;
#     if ($yesno[$i] == 1){$tsum++;}
#     if ($tsum > 0){
#       print "$scores[$i]\t",($yesrateabove[scalar(@sorted)-1-$n])/($tsum/$n)-1,"\n";
#     }
#   }

  
}


if ($answersscatterplot eq'yes' or $truerate eq 'yes'){
  my (%yesrate,%trueyesrate,%rate);
  my (@yrate,@trate);
open F,("$answers");
while (<F>){
  unless (/^\S*\s+ont/){
    @temp=split /\t/,$_;chomp($temp[scalar(@temp)-1]);
    if ($temp[8] =~ /\S+/ or 1 == 1){
    $rate{$temp[1]}++;
    if ($temp[4] eq 'True'){
      $yesrate{$temp[1]}++;
      if ($temp[8] eq 'True'){
      $trueyesrate{$temp[1]}++;
    }
    }
  }
  }
}
  close F;
  foreach my $i (keys(%yesrate)){
    unless (exists($trueyesrate{$i})){$trueyesrate{$i}=0;}
            $trueyesrate{$i}=$trueyesrate{$i}/$yesrate{$i};
    $yesrate{$i}=$yesrate{$i}/$rate{$i};
    unless ($truerate eq 'yes'){
      print "$yesrate{$i}\t$trueyesrate{$i}\n";
    }
    push @yrate,$yesrate{$i};
    push @trate,$trueyesrate{$i};
  }

  my $n=0;my $tsum=0;
  my @sorted = OrderArray(@yrate);
  foreach my $i (@sorted){
    $n++;
    $tsum+=$trate[$i];
    my $averate=$tsum/$n;
    print "$yrate[$i]\t$averate\n";
  }
  
}


if ($scoretruerate eq 'yes'){
  my (@trues,@scores);
open F,("$answers");
while (<F>){
  unless (/^\S*\s+ont/){
    @temp=split /\t/,$_;chomp($temp[scalar(@temp)-1]);
    unless (exists($termyesrates{$temp[1]})){
$termyesrates{$temp[1]}=0;
    }
    unless ($termyesrates{$temp[1]} >= $ratefilter){
    if ($temp[8] =~ /\S+/ or 1 == 1 ){
      push @scores,$temp[5];
      if ($temp[4] eq 'True'){
push @trues,1;
      }
      else{
push @trues,0;
      }
    }
  }
  }
}
  close F;
  my @sorted = OrderArray(@scores);
  my $i=0;my $true=0;my $allyes=0;my $n=0;
  foreach my $j (0 .. scalar(@sorted)-1){
    $j = $sorted[scalar(@sorted)-1-$j];
    $i++;
    $n++;
    if ($trues[$j] == 1){
      $allyes++;
$true++;
}
    if ($i == $scorebinsize){
      my $rate=$true/$i;
      my $cumrate=$allyes/$n;
#      print "$n\t$cumrate\n";
#      print "$scores[$j]\t$cumrate\n";
      print "$scores[$j]\t$n\n";
            $true=0;$i=0;
    }
  }
  my $yr=$allyes/$n;
#print "\n0\t$yr\n$n\t$yr\n";
print "\n0\t$yr\n1\t$yr\n";
}


if ($decoyscatter eq 'yes'){
  my (%yes,%no,%decoy,%pred);
open F,("$answers");
while (<F>){
  unless (/^\S*\s+ont/){
    @temp=split /\t/,$_;chomp($temp[scalar(@temp)-1]);
    if ($temp[4] eq 'True'){
$yes{$temp[1]}++;
    }
    elsif ($temp[4] eq 'False'){
$no{$temp[1]}++;
    }
    if ($temp[8] eq 'True' and $temp[4] eq 'True'){
$pred{$temp[1]}++;
    }
    elsif ($temp[8] eq 'False' and $temp[4] eq 'True'){
$decoy{$temp[1]}++;
    }
  }
}
  close F;
  my (@yesrate,@decoyrate);
  foreach my $term (keys(%yes)){
    unless (exists($yes{$term})){$yes{$term} = 0}    unless (exists($decoy{$term})){$decoy{$term} = 0}
    unless (exists($no{$term})){$no{$term} = 0}    unless (exists($pred{$term})){$pred{$term} = 0}
#    print $yes{$term}/($yes{$term}+$no{$term}),"\t",$decoy{$term}/($decoy{$term}+$pred{$term}),"\n";
    print $yes{$term}+$no{$term},"\t",$decoy{$term},"\n";
#    print $yes{$term}/($yes{$term}+$no{$term}),"\t",$pred{$term},"\n";
#    push @decoyrate,$decoy{$term}/($decoy{$term}+$pred{$term});
    push @decoyrate,$decoy{$term};
#    if ($decoy{$term} > 250){print STDERR "$term\t$decoy{$term}\n";}
#    push @yesrate,$yes{$term}/($yes{$term}+$no{$term});
    push @yesrate,$yes{$term}+$no{$term};
  }
  print "\n";
  my $n=0;my $tsum=0;
  my @sorted = OrderArray(@yesrate);
  foreach my $i (0 .. scalar(@sorted)-1){
    $i = $sorted[scalar(@sorted)-1-$i];
#    $i = $sorted[$i];
    $n++;
    $tsum+=$decoyrate[$i];
    my $averate=$tsum/$n;
    print "$yesrate[$i]\t$averate\n";
  }
  
}





if ($permutation_byquestion_top_5 eq 'yes'){
  my $type = 'sum';
  my (%term,%isyes,%score,%predicted,%hasyes);
open F,("$answers");
while (<F>){
  unless (/^\S*\s+ont/){
    @temp=split /\t/,$_;chomp($temp[scalar(@temp)-1]);
    if ($temp[8] =~ /\S+/ or 1 == 1){#include recall or not
            unless (exists($termyesrates{$temp[1]})){$termyesrates{$temp[1]}=0;}
            unless (exists($termyesrates_norecall{$temp[1]})){$termyesrates_norecall{$temp[1]}=0;}
      if (exists($term{$temp[1]})){$term{$temp[1]}=$term{$temp[1]}.",$temp[2]";}else{$term{$temp[1]}=$temp[2]}
if ($temp[5] >= 0){	$score{$temp[1]}{$temp[2]}=$temp[5];}else{$score{$temp[1]}{$temp[2]}=0;}
    if ($temp[4] eq 'True'){$isyes{$temp[1]}{$temp[2]}=1;$hasyes{$temp[1]}=1;}elsif($temp[4] eq 'False'){$isyes{$temp[1]}{$temp[2]}=0;}
    }
  }
}
  close F;

  for my $rate (0 .. 500){
  my @obstop;
    my ($mean,$variance,$zscore,$n);
    my $sums=0;
  my $average=0;
    $top=0.01+0.001*$rate;
    $ratefilter=0.05;
    my %permutes;    my %totalsums; my %binaries;my %binary;
  #generate permutation data
  foreach my $t (keys(%term)){
    if ( $termyesrates{$t} < $ratefilter){
    my @temp=split /,/,$term{$t};
    my @a;my @s;
    foreach my $j (@temp){
       push @a,$isyes{$t}{$j};
      push @s,$score{$t}{$j};
    }
    my ($r,$i,$sum,$bin);my @rand;
    for my $j (1 .. $perms){
       @rand = shuffle(@a);
      $i=0;$sum=0;$bin=0;
      foreach $r (@rand){
	if ($r == 1){$sum+=$s[$i]};
	if ($r == 1 and $s[$i] >= $threshold){$bin++;}
	$i++;
      }
      $sums+=$sum;
      unless (exists($permutes{$j}{$t})){$permutes{$j}{$t}=0;}
      unless (exists($binaries{$j}{$t})){$binaries{$j}{$t}=0;}
       $permutes{$j}{$t}+=$sum;
       $totalsums{$j}+=$sum;
	 $binary{$j}+=$bin;
       $binaries{$j}{$t}+=$bin;
     }
  }
  }

    $average=$sums;
    

    my @obsscores;
  
  #calculate observed
  my $tot=0;my $tms=0;
   foreach my $t (keys(%term)){
       if ( $termyesrates{$t} < $ratefilter){
	 $tms++;
     my @temp=split /,/,$term{$t};
     my $tt=0;
        foreach my $j (@temp){
	  if ($isyes{$t}{$j} == 1){
	    if ($type eq 'sum'){
	      $tt+=$score{$t}{$j};
	    }
	      else{
		if ($score{$t}{$j} >= $threshold){$tt++;}
	      }
	  }
	  	    	    $tot+=$tt;
	}
     push @obsscores,$tt;
   }
   }


      #calculate per permute
  for my $j (1 .. $perms){
    my @p;
    foreach my $t (keys(%term)){
    if ( $termyesrates{$t} < $ratefilter){
	    if ($type eq 'sum'){
	      push @p,$permutes{$j}{$t};
	    }
	    else{
	      push @p,$binaries{$j}{$t};
	    }
  }
    }
    my @sorted = OrderArray(@p);
my $top5 = ($p[$sorted[int(scalar(@sorted)*(1-$top)+0.5)-1]]+$p[$sorted[int(scalar(@sorted)*(1-$top)+0.5)-1+1]])/2;
push @obstop,$top5;
  }


    my @sorted = OrderArray(@obsscores);
  my %dist;my $min=scalar(@obsscores);my $max=0;
    foreach my $stop (@obstop){
      my $nn=0;
    foreach my $sc (0 .. scalar(@sorted)-1){
      $sc=$obsscores[scalar(@sorted)-1-$sc];
	if ($sc > $stop){
	  $nn++;
	}
    }
$dist{$nn}++;
	  if ($nn > $max){$max = $nn;}elsif($nn < $min){$min=$nn;}
      }
    

    my $nnn=0;my $ttt=0;my $vvv=0;
    for my $i ($min .. $max){
      unless (exists($dist{$i})){$dist{$i}=0;}
      $nnn+=$dist{$i};
      $vvv+=$dist{$i}*$i**2;
      $ttt+=$i*$dist{$i};
print "$i\t$dist{$i}\n";
    }
    print "\n";
    
#    print $tms*$top,"\t$perms\t$tms\n";

 
       $mean=$ttt/$nnn;
       $variance=$vvv/$nnn-$mean**2;
       $variance=$variance**0.5;
     	$zscore=($mean-$tms*$top)/$variance;
  print  100*(($ttt/$nnn)/($tms*$top)-1),"\%  correct: ",$ttt/$nnn-$tms*$top,"     total: ",$ttt/$nnn,"   top: $top    Zscore:  $zscore\n";


}
  
}





if ($recallyesrates eq 'yes'){
  my (%rqs,%rps);
open F,("$answers");
while (<F>){
  unless (/^\S*\s+ont/){
    @temp=split /\t/,$_;chomp($temp[scalar(@temp)-1]);
    unless ($temp[8] =~ /\S+/){
      $rqs{$temp[1]}=1;
      $rps{$temp[2]}=1;
    }
  }
}
  close F;

my $rpy=0;my $rpn=0;my $rqy=0;my $rqn=0;my $rrpy=0;my $rrpn=0;my $rrqy=0;my $rrqn=0;my $yes=0;my $no=0;
open F,("$answers");
while (<F>){
  unless (/^\S*\s+ont/){
    @temp=split /\t/,$_;chomp($temp[scalar(@temp)-1]);
    if ($temp[8] =~ /\S+/){
      if ($temp[4] eq 'True'){
	$yes++;
			      if (exists($rps{$temp[2]})){$rpy++;}
  			      if (exists($rqs{$temp[1]})){$rqy++;}
    }
      else{
	$no++;
			      if (exists($rps{$temp[2]})){$rpn++;}
  			      if (exists($rqs{$temp[1]})){$rqn++;}
      }
    }
    else{
   if ($temp[4] eq 'True'){
			      if (exists($rps{$temp[2]})){$rrpy++;}else{die;}
  			      if (exists($rqs{$temp[1]})){$rrqy++;}else{die;}

      }
      else{
			      if (exists($rps{$temp[2]})){$rrpn++;}else{die;}
  			      if (exists($rqs{$temp[1]})){$rrqn++;}else{die;}
      }
    }
  }
}
  close F;
  print "Yes rate before recall: ",$yes/($yes+$no),"\n";
print "Yes rate for recall people first time: ",$rpy/($rpy+$rpn),"\nYes rate for people in recall second time: ",$rrpy/($rrpy+$rrpn),"\n";
print "Yes rate for recall questions first time: ",$rqy/($rqy+$rqn),"\nYes rate for recall questions second time: ",$rrqy/($rrqy+$rrqn),"\n";
}


if ($permutation_test_question eq 'yes'){
  my (%term,%isyes,%score,%predicted);
open F,("$answers");
while (<F>){
  unless (/^\S*\s+ont/){
     @temp=split /\t/,$_;chomp($temp[scalar(@temp)-1]);
    if ($temp[8] =~ /\S+/ or 1 == 1){#include recall or not
           unless (exists($termyesrates{$temp[1]})){$termyesrates{$temp[1]}=0;}
            unless (exists($termyesrates_norecall{$temp[1]})){$termyesrates_norecall{$temp[1]}=0;}
      if (exists($term{$temp[1]})){$term{$temp[1]}=$term{$temp[1]}.",$temp[2]";}else{$term{$temp[1]}=$temp[2]}
if ($temp[5] >= 0){	$score{$temp[1]}{$temp[2]}=$temp[5];}else{$score{$temp[1]}{$temp[2]}=0;}
    if ($temp[4] eq 'True'){$isyes{$temp[1]}{$temp[2]}=1;}elsif($temp[4] eq 'False'){$isyes{$temp[1]}{$temp[2]}=0;}
    }
  }
}
  close F;

  for my $rate (1 .. 100){
    my ($mean,$variance,$zscore,$n);
    my $sums=0;
    my $average=0;
    $ratefilter=$rate/1000;
    $threshold=0+$rate/1000;
    $ratefilter=0.05;
    $threshold=0.023;
    my %permutes;    my %totalsums; my %above; my %aboveyes;
  #generate permutation data
  foreach my $t (keys(%term)){
    if ( $termyesrates{$t} < $ratefilter){
    my @temp=split /,/,$term{$t};
    my @a;my @s;
    my $testrunning=0;
    foreach my $j (@temp){
       push @a,$isyes{$t}{$j};
      push @s,$score{$t}{$j};
    }
    my ($r,$i,$sum);$permutes{$t}='';my @rand;
           my $testsa=0;my $testsb=0;my $testa=0;my $testb=0;
    for my $j (1 .. $perms){
       @rand = shuffle(@a);
       $i=0;$sum=0;
       foreach $r (@rand){
	if ($s[$i] >= $threshold){$above{$j}++;$testa++;}else{$testb++;}
	if ($r == 1){
	  $sum+=$s[$i];
	  if ($s[$i] >= $threshold){$aboveyes{$j}++;$testsa++;}else{$testsb++;}
	}
	$i++;
      }
       $sums+=$sum;
       $permutes{$t}.=",$sum";
       $totalsums{$j}+=$sum;
     }
    unless ($testa == 0 or $testb == 0){
#    print $testsa/$testa-$testsb/$testb,"\t",$testsa/$testa,"\t",$testsb/$testb,"\t$testsa $testa $testsb $testb\t$testrunning\n";
    $testrunning=$testrunning+$testsa/$testa-$testsb/$testb;
  }
  }
  }
    $average=$sums;

  
  #calculate observed
      my $yesses=0;
    my $answers=0;
    my $aboveobs=0;
    my $aboveyesobs=0;
my $tot=0; 
   foreach my $t (keys(%term)){
   if ( $termyesrates{$t} < $ratefilter){
     my @temp=split /,/,$term{$t};
     foreach my $j (@temp){
       $answers++;
  	if ($score{$t}{$j} >= $threshold){$aboveobs++;}
     if ($isyes{$t}{$j} == 1){
  	if ($score{$t}{$j} >= $threshold){$aboveyesobs++;}
	 $yesses++;
      $tot+=$score{$t}{$j};
    }
   }
 }
   }
    $average=$average/$perms;
    my $ratio;
	    my $obsrateabove=$aboveyesobs/$aboveobs;
	    my $obsratebelow=($yesses-$aboveyesobs)/($answers-$aboveobs);
    my $obscorrect=$obsrateabove-$obsratebelow;

    
    if (1 == 0){#this is for sum of scores
      $n=0;
      for my $j (1 .. $perms){
	$mean+=$totalsums{$j};
	$variance+=$totalsums{$j}**2;
	$n++;
#print "$totalsums{$j}\n";
      }
      $mean=$mean/$n;
      $variance=$variance/$n-$mean**2;
      $variance=$variance**0.5;
      $zscore=($tot-$mean)/$variance;
    print  "$ratefilter\t$zscore\t$ratio\t$tot\t$average\t$mean\t$variance\n";
   }

    if (1 == 1){#this is for yes rate above threshold
      $n=0;
      for my $j (1 .. $perms){
	unless (exists($above{$j})){$above{$j}=0;}unless (exists($aboveyes{$j})){$aboveyes{$j}=0;}
	my $rateabove=$aboveyes{$j}/$above{$j};
	    my $ratebelow=($yesses-$aboveyes{$j})/($answers-$above{$j});
	    my $correct=$rateabove-$ratebelow;
	$mean+=$correct;
	$variance+=$correct**2;
	$n++;
#print "$correct\t$aboveyes{$j}\t",$yesses-$aboveyes{$j},"\t$rateabove\t$ratebelow\t$obsrateabove\t$obsratebelow\n";
      }
      $mean=$mean/$n;
      $variance=$variance/$n-$mean**2;
      $variance=$variance**0.5;
      $zscore=($obscorrect-$mean)/$variance;
    print  "$threshold\t$ratefilter\t$zscore\t$obscorrect\t$mean\t$variance\n";
  }


    
}
  
}


if ($permutation_test_person eq 'yes'){
  my (%term,%isyes,%score,%predicted);
open F,("$answers");
while (<F>){
  unless (/^\S*\s+ont/){
        @temp=split /\t/,$_;chomp($temp[scalar(@temp)-1]);
      if ($temp[8] =~ /\S+/ or 1 == 1){#include recall or not
           unless (exists($termyesrates{$temp[1]})){$termyesrates{$temp[1]}=0;}
            unless (exists($termyesrates_norecall{$temp[1]})){$termyesrates_norecall{$temp[1]}=0;}
      #to reverse people and questions:
      my $te=$temp[1];$temp[1]=$temp[2];$temp[2]=$te;
      if (exists($term{$temp[1]})){$term{$temp[1]}=$term{$temp[1]}.",$temp[2]";}else{$term{$temp[1]}=$temp[2]}
if ($temp[5] >= 0){	$score{$temp[1]}{$temp[2]}=$temp[5];}else{$score{$temp[1]}{$temp[2]}=0;}
    if ($temp[4] eq 'True'){$isyes{$temp[1]}{$temp[2]}=1;}elsif($temp[4] eq 'False'){$isyes{$temp[1]}{$temp[2]}=0;}
    }
  }
}
  close F;

  for my $rate (10 .. 1000){
    my ($mean,$variance,$zscore,$n);
    my $sums=0;
  my $average=0;
$ratefilter=$rate/1000;
    my %permutes;    my %totalsums;
  #generate permutation data
  foreach my $t (keys(%term)){
    my @temp=split /,/,$term{$t};
    my @a;my @s;
    foreach my $j (@temp){
    if ( $termyesrates{$t} < $ratefilter){
       push @a,$isyes{$t}{$j};
       push @s,$score{$t}{$j};
  }
    }
    my ($r,$i,$sum);$permutes{$t}='';my @rand;
    for my $j (1 .. $perms){
       @rand = shuffle(@a);
      $i=0;$sum=0;
      foreach $r (@rand){
	if ($r == 1){$sum+=$s[$i]};
	$i++;
      }
       $sums+=$sum;
       $permutes{$t}.=",$sum";
       $totalsums{$j}+=$sum;
     }
  }

    $average=$sums;

  #calculate observed
  my $tot=0;
   foreach my $t (keys(%term)){
     my @temp=split /,/,$term{$t};
        foreach my $j (@temp){
   if ( $termyesrates{$t} < $ratefilter){
	  if ($isyes{$t}{$j} == 1){
	    $tot+=$score{$t}{$j};
	  }
	}
	}
   }
    $average=$average/$perms;
        my $ratio=$tot/$average;

    if (1 == 1){
      for my $j (1 .. $perms){
	$mean+=$totalsums{$j};
	$variance+=$totalsums{$j}**2;
	$n++;
#print "$totalsums{$j}\n";
      }
      $mean=$mean/$n;
      $variance=$variance/$n-$mean**2;
      $variance=$variance**0.5;
      	$zscore=($tot-$mean)/$variance;
    }

    print  "$ratefilter\t$zscore\t$ratio\t$tot\t$average\t$mean\t$variance\n";

}
  
}





if ($permutation_test_decoyszero eq 'yes'){
  my (%term,%isyes,%score,%predicted,%decoy);my (@allscores,@terms);
    my $decoyyes=0;my $decoyans=0;
open F,("$answers");
while (<F>){
  unless (/^\S*\s+ont/){
    @temp=split /\t/,$_;chomp($temp[scalar(@temp)-1]);
    if ($temp[5] >= $zero){#exclude decoy set of low scores
       if ($temp[8] eq 'True'){$decoy{$temp[1]}=0;if ($temp[4] eq 'True'){$decoy{$temp[1]}=1;}}
          unless (exists($termyesrates{$temp[1]})){$termyesrates{$temp[1]}=0;}
            unless (exists($termyesrates_norecall{$temp[1]})){$termyesrates_norecall{$temp[1]}=0;}
      if (exists($term{$temp[1]})){$term{$temp[1]}=$term{$temp[1]}.",$temp[2]";}else{$term{$temp[1]}=$temp[2];}
if ($temp[5] >= 0){push @allscores,$temp[5];$score{$temp[1]}{$temp[2]}=$temp[5];}else{push @allscores,0;$score{$temp[1]}{$temp[2]}=0;}
    if ($temp[4] eq 'True'){$isyes{$temp[1]}{$temp[2]}=1;}elsif($temp[4] eq 'False'){$isyes{$temp[1]}{$temp[2]}=0;}
     }
    else{
      if ($temp[4] eq 'True'){
$decoyyes++;
      }
$decoyans++;
       }
  }
}
  close F;
  foreach my $t (keys(%term)){
      push @terms,$t;
  }

      my $testrunning=0;my $nn=0;
  for my $rate ($arg1 .. $arg2){
    my $sums=0;
    my $average=0;
    $ratefilter=$rate/100;
    $threshold=0.018+$rate/1000;
#    $ratefilter=0.05;
    $threshold=0.023;
 my %totalsums; my %above; my %aboveyes;
my @randscores;
    for my $j (1 .. $perms){
           my $testsa=0;my $testsb=0;my $testa=0;my $testb=0; my $sum=0;
    #generate permutation data
my @randscores=shuffle(@allscores);
	   foreach my $t (@terms){
   if ( $termyesrates{$t} < $ratefilter){
    my @person=split /,/,$term{$t};
       foreach my $p (@person){
my	 $r=pop(@randscores);
	 my $a=$isyes{$t}{$p};
	 	if ($r >= $threshold){$above{$j}++;$testa++;}else{$testb++;}
	if ($a == 1){
	  $sum+=$r;
	  if ($r >= $threshold){$aboveyes{$j}++;$testsa++;}else{$testsb++;}
	}
       }
           $sums+=$sum;
  }
 }
       $totalsums{$j}=$sum;
	   if ($testb > 0 and $testa > 0){
	     $testrunning+=$testsa/$testa-$testsb/$testb;$nn++;
	   }
    }

    
  #calculate observed
      my $yesses=0;
    my $answers=0;
    my $aboveobs=0;
    my $aboveyesobs=0;
    my $tot=0;
    my %include;
   foreach my $t (keys(%term)){
     if ( $termyesrates{$t} < $ratefilter){
      $include{$t}=1;
     my @temp=split /,/,$term{$t};
     foreach my $j (@temp){
	 $answers++;
  	if ($score{$t}{$j} >= $threshold){$aboveobs++;}
     if ($isyes{$t}{$j} == 1){
  	if ($score{$t}{$j} >= $threshold){$aboveyesobs++;}
	 $yesses++;
      $tot+=$score{$t}{$j};
    }
   }
 }
   }
    $average=$sums/$perms;
        my $ratio=$tot/$average;
	  
    
    if (1 == 1){#this is for sum of scores
     my ($mean,$variance,$zscore,$n);
     $n=0;
      for my $j (1 .. $perms){
	$mean+=$totalsums{$j};
	$variance+=$totalsums{$j}**2;
	$n++;
#print "$totalsums{$j}\n";
      }
      $mean=$mean/$n;
      $variance=$variance/$n-$mean**2;
      $variance=$variance**0.5;
      $zscore=($tot-$mean)/$variance;
    print  "$ratefilter\t$zscore\t$ratio\t$tot\t$average\t$mean\t$variance\n";
   }

    if (1 == 1){#this is for yes rate above threshold
     my ($mean,$variance,$zscore,$n);
   my $obsrateabove=$aboveyesobs/$aboveobs;
#	    my $obsratebelow=($yesses-$aboveyesobs)/($answers-$aboveobs);
my $obsratebelow=$decoyyes/$decoyans;
    my $obscorrect=$obsrateabove-$obsratebelow;
    $n=0;
      for my $j (1 .. $perms){
	unless (exists($above{$j})){$above{$j}=0;}unless (exists($aboveyes{$j})){$aboveyes{$j}=0;}
	my $rateabove=$aboveyes{$j}/$above{$j};
#	my $ratebelow=($yesses-$aboveyes{$j})/($answers-$above{$j});
	my $ratebelow=$decoyyes/$decoyans;
	    my $correct=$rateabove-$ratebelow;
	$mean+=$correct;
	$variance+=$correct**2;
	$n++;
print "$correct\t$aboveyes{$j}\t",$yesses-$aboveyes{$j},"\t$rateabove\t$ratebelow\t$obsrateabove\t$obsratebelow\n";
      }
      $mean=$mean/$n;
      $variance=$variance/$n-$mean**2;
      $variance=$variance**0.5;
     $zscore=($obscorrect-$mean)/$variance;
    print  "$threshold\t$ratefilter\t$zscore\t$obscorrect\t$mean\t$variance\t",$obsrateabove/$obsratebelow,"\t$obsrateabove\t",($yesses-$aboveyesobs)/($answers-$aboveobs),"\t",scalar(keys(%include))/scalar(keys(%term)),"\t",$decoyyes/$decoyans,"\n";
  }
  }
}



if ($permutation_test_global eq 'yes'){
  my (%term,%isyes,%score,%predicted,%decoy);my (@allscores,@terms);
open F,("$answers");
while (<F>){
  unless (/^\S*\s+ont/){
    @temp=split /\t/,$_;chomp($temp[scalar(@temp)-1]);
    if ($temp[8] =~ /\S+/ or 1 == 1){#include recall or not
       if ($temp[8] eq 'True'){$decoy{$temp[1]}=0;if ($temp[4] eq 'True'){$decoy{$temp[1]}=1;}}
          unless (exists($termyesrates{$temp[1]})){$termyesrates{$temp[1]}=0;}
            unless (exists($termyesrates_norecall{$temp[1]})){$termyesrates_norecall{$temp[1]}=0;}
      if (exists($term{$temp[1]})){$term{$temp[1]}=$term{$temp[1]}.",$temp[2]";}else{$term{$temp[1]}=$temp[2];}
if ($temp[5] >= 0){push @allscores,$temp[5];$score{$temp[1]}{$temp[2]}=$temp[5];}else{push @allscores,0;$score{$temp[1]}{$temp[2]}=0;}
    if ($temp[4] eq 'True'){$isyes{$temp[1]}{$temp[2]}=1;}elsif($temp[4] eq 'False'){$isyes{$temp[1]}{$temp[2]}=0;}
    }
  }
}
  close F;
  foreach my $t (keys(%term)){
      push @terms,$t;
  }
if (defined($ARGV[1])){$arg1=$ARGV[1];}else{$arg1=1;}
if (defined($ARGV[2])){$arg2=$ARGV[2];}else{$arg2=1;}

      my $testrunning=0;my $nn=0;
  for my $rate ($arg1 .. $arg2){
    my $decoyyes=0;my $decoyans=0;
    my $sums=0;
    my $average=0;
    $ratefilter=$rate/100;
    $threshold=0.021+$rate/10000;
#    $ratefilter=0.05;
    $threshold=0.022;
 my %totalsums; my %above; my %aboveyes;
my @randscores;
    for my $j (1 .. $perms){
           my $testsa=0;my $testsb=0;my $testa=0;my $testb=0; my $sum=0;
    #generate permutation data
my @randscores=shuffle(@allscores);
	   foreach my $t (@terms){
   if ( $termyesrates{$t} < $ratefilter){
    my @person=split /,/,$term{$t};
       foreach my $p (@person){
my	 $r=pop(@randscores);
	 my $a=$isyes{$t}{$p};
	 	if ($r >= $threshold){$above{$j}++;$testa++;}else{$testb++;}
	if ($a == 1){
	  $sum+=$r;
	  if ($r >= $threshold){$aboveyes{$j}++;$testsa++;}else{$testsb++;}
	}
       }
           $sums+=$sum;
  }
 }
       $totalsums{$j}=$sum;
	   if ($testb > 0 and $testa > 0){
	     $testrunning+=$testsa/$testa-$testsb/$testb;$nn++;
	   }
    }

    
  #calculate observed
      my $yesses=0;
    my $answers=0;
    my $aboveobs=0;
    my $aboveyesobs=0;
    my $tot=0;
    my %include;
   foreach my $t (keys(%term)){
     if ( $termyesrates{$t} < $ratefilter){
      $include{$t}=1;
       if (defined($decoy{$t})){
	 $decoyans++;
	 if ($decoy{$t}==1){$decoyyes++;}
       }
     my @temp=split /,/,$term{$t};
     foreach my $j (@temp){
	 $answers++;
  	if ($score{$t}{$j} >= $threshold){$aboveobs++;}
     if ($isyes{$t}{$j} == 1){
  	if ($score{$t}{$j} >= $threshold){$aboveyesobs++;}
	 $yesses++;
      $tot+=$score{$t}{$j};
    }
   }
 }
   }
    $average=$sums/$perms;
        my $ratio;
unless ($average != 0){$ratio=0;}else{$ratio=$tot/$average;}

    
    if (1 == 1){#this is for sum of scores
     my ($mean,$variance,$zscore,$n);
     $n=0;
      for my $j (1 .. $perms){
	$mean+=$totalsums{$j};
	$variance+=$totalsums{$j}**2;
	$n++;
#print "$totalsums{$j}\n";
      }
      $mean=$mean/$n;
      $variance=$variance/$n-$mean**2;
      $variance=$variance**0.5;
if ($tot-$mean == 0){$zscore=0;}else{      $zscore=($tot-$mean)/$variance;}
    print  "$ratefilter\t$zscore\t$ratio\t$tot\t$average\t$mean\t$variance\n";
   }

    if (1 == 1){#this is for yes rate above threshold
     my ($mean,$variance,$zscore,$n);
   my $obsrateabove=$aboveyesobs/$aboveobs;
#	    my $obsratebelow=($yesses-$aboveyesobs)/($answers-$aboveobs);
my $obsratebelow=$decoyyes/$decoyans;
     #    my $obscorrect=$obsrateabove-$obsratebelow;
     my $obscorrect=$aboveyesobs;
    $n=0;
      for my $j (1 .. $perms){
	unless (exists($above{$j})){$above{$j}=0;}unless (exists($aboveyes{$j})){$aboveyes{$j}=0;}
	my $rateabove=$aboveyes{$j}/$above{$j};
#	my $ratebelow=($yesses-$aboveyes{$j})/($answers-$above{$j});
	my $ratebelow=$decoyyes/$decoyans;
	#	    my $correct=$rateabove-$ratebelow;
	my $correct=$aboveyes{$j};
	$mean+=$correct;
	$variance+=$correct**2;
	$n++;
print "$totalsums{$j}\t$correct\t$aboveyes{$j}\t",$yesses-$aboveyes{$j},"\t$rateabove\t$ratebelow\t$obsrateabove\t$obsratebelow\n";
      }
      $mean=$mean/$n;
      $variance=$variance/$n-$mean**2;
      $variance=$variance**0.5;
     $zscore=($obscorrect-$mean)/$variance;
    print  "$threshold\t$ratefilter\t$zscore\t$obscorrect\t$mean\t$variance\t",$obsrateabove/$obsratebelow,"\t$obsrateabove\t",($yesses-$aboveyesobs)/($answers-$aboveobs),"\t",scalar(keys(%include))/scalar(keys(%term)),"\t",$decoyyes/$decoyans,"\n";
  }


    
} 
}



#SUB-ROUTINE-----------------------------------------------
sub OrderArray{
#  This reads in an array of values and orders them 
#returning the list of the order starting from lowest going to highest

my @values=@_;
my %map;
for my $i (0 .. scalar(@values)-1){
$map{$i}=$values[$i];
}
my @listout = sort { $map{$a} <=> $map{$b} } keys %map; 

return @listout;
}
#----------------------------------------------------------

