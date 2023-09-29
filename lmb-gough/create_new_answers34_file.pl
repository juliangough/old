#!/usr/bin/perl -w
use strict;


my $answers='/data2/gough/clu/answers34.tsv';
my $datadir='/home/snow/data3/impute';
my %scores;

for my $i (0 .. 4){
  my %mapping;
  open M,("/home/snow/data3/impute/data$i/temp/predictor.input.vcf");
  while (<M>){
    if (/#CHROM\s+POS.+INFO\s+FORMAT\s+(\S+.+\S+)$/){
my @list=split /\s+/,$1;
for my $j (1 .. scalar(@list)){
  $mapping{$j}=$list[$j-1];
}
last;
}
  }
  close M;
  open D,("$datadir/data$i/temp/predictor.allscores.tsv");
  while (<D>){
    my @temp=split /\s+/,$_;chomp($temp[scalar(@temp)-1]);
    if ($temp[1] > 2504){
      my $id=$mapping{$temp[1]};
    if (exists( $scores{$temp[0]}{$temp[1]})){
print STDERR "WARNING: score repeated for $temp[0] $id in $i old:$scores{$temp[0]}{$id} new:$temp[3]:$temp[6]\n";
    }
    else{
      $scores{$temp[0]}{$id}="$temp[3]:$temp[6]";
    }
  }
}
  close D;
#   open LS,("ls $datadir/data$i/temp/*.scores |");
#   while (<LS>){
#     if (/([^\/]+)\.scores$/){
#       my $term=$1;
#       open D,("$_");
#       while (<D>){
# 	unless (/^#/){
# 	  my @temp=split /\s+/,$_;chomp($temp[scalar(@temp)-1]);
# 	  if ($temp[0] > 2504){
# 	my $id=$mapping{$temp[0]};
# 	if (exists($scores{$term}{$id})){
# 	  $scores{$term}{$id} =~ /(\S+):(\S+)/;
# 	  unless ($2 == $temp[8]){
# print STDERR "WARNING different scores in different files for $id ($temp[0]) for term $term where $term.scores has $temp[8] and predictor.allscores.tsv has $1\n";
# }
# 	  else{
# #print STDERR "GOOD same scores in different files for $id ($temp[0]) for term $term where $term.scores has $temp[8] and predictor.allscores.tsv has $1\n";
# 	  }
# 	}
# 	else{
# print STDERR "WARNING score missing for $id ($temp[0]) for term $term where $term.scores has $temp[8] and predictor.allscores.tsv has nothing\n";
# 	}
#       }
#       }
#       }
#       close D;
#     }
#   }
#   close LS;
  #   open LS,("ls $datadir/data$i/*.predictions.tsv |");
#   while (<LS>){
#     if (/(\d+)\.predictions\.tsv$/){
#       my $id=$1;
#       open D,(" $datadir/data$i/$id.predictions.tsv");
#       while (<D>){
# 	my @temp=split /\t/,$_;
# 	if (exists($scores{$temp[2]}{$id})){
# 	  $scores{$temp[2]}{$id} =~ /(\S+):(\S+)/;
# 	  unless ($1 == $temp[1]){
# print STDERR "WARNING different scores in different files for $id in $i   file:$id.predictor.tsv has $temp[1] and predictor.allscores.tsv has $1\n";
# }
# 	  else{
# print STDERR "GOOD same scores in different files for $id in $i   file:$id.predictor.tsv has $temp[1] and predictor.allscores.tsv has $1\n";
# 	  }
# 	}
# 	else{
# 	  $scores{$temp[2]}{$id}="$temp[1]:0";
# 	}
#       }
#       close D;
#     }
#   }
#   close LS;
}

open F,("$answers");
while (<F>){
  unless (/ont\s+uniqueid/){
    my @temp=split /\s+/,$_;chomp($temp[scalar(@temp)-1]);
    if (exists( $scores{$temp[1]}{$temp[2]})){
      if ($scores{$temp[1]}{$temp[2]} =~ /^(\S+):(\S+)$/){
	print "$temp[0]\t$temp[1]\t$temp[2]\t$temp[3]\t$temp[4]\t$1\t$2\t$temp[7]\t$temp[8]\n";
	if (($temp[5]-$1)**2 > 0.022**2){print STDERR "Big change for $temp[1] $temp[2] old:$temp[5] new:$1\n";}
      }
      else{die;}
    }
    else{
      print STDERR "WARNING: score missing for $temp[1] $temp[2] so setting to zero\n";
	print "$temp[0]\t$temp[1]\t$temp[2]\t$temp[3]\t$temp[4]\t0.0\t0.0\t$temp[7]\t$temp[8]\n";
	if (($temp[5])**2 > 0.022**2){print STDERR "Big change for $temp[1] $temp[2] old:$temp[5] new:missing - set to zero\n";}
    }
  }
}
close F;
