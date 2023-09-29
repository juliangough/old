#!/usr/bin/perl -w
use strict;

my $dir;
if (defined($ARGV[0])){$dir=$ARGV[0]}else{print STDERR "Usage: remove_vcf_duplicates.pl [dir]\n";die;}

open LS,("ls $dir |");
while (<LS>){
  my $file = $_;chomp($file);
  my @lines;my %line;my %score;
  open F,("$file");
  while (<F>){
    if(/^(\S+)\t(\d+)\t(\S+)\t\S+\t\S+\t\S+\t\S+\t\S+\t\S+\t(\S.+\S)$/){
      my $loc="$1:$2";my $new=$_;my $rs=$3;
      push @lines,$loc;
      my @temp=split /\t/,$4;my $x=0;my $y=0;
      foreach my $a (@temp){	if ($a =~ /[01]/){$x++;}	else{$y++;	}      }
      if (exists($line{$loc})){
	if ($line{$loc} =~ /^(\S+)\t(\d+)\t(\S+)\t\S+\t\S+\t\S+\t\S+\t\S+\t\S+\t(\S.+\S)$/){
      my @temp2=split /\t/,$4;my $x2=0;my $y2=0;my $rs2=$3;
      foreach my $a (@temp2){	if ($a =~ /[01]/){$x2++;}	else{$y2++;	}      }
      if (($x2-$y2) == $score{$loc}){
	if ($rs =~ /^i/ and $rs2 =~ /^rs/){
	print STDERR "Replace($rs): $line{$loc}New($rs2):     $new";
	$line{$loc}=$new;
	 $score{$loc}=$x2-$y2;
	}
	else{
	print STDERR "Keep($rs):    $line{$loc}Reject($rs2): $new";
	}
      }
      elsif(($x2-$y2) > $score{$loc}){
	print STDERR "Replace($score{$loc}): $line{$loc}New(",($x2-$y2),"):     $new";
	$line{$loc}=$new;
	 $score{$loc}=$x2-$y2;
      }
      else{
	print STDERR "Keep($score{$loc}):    $line{$loc}Reject(",($x2-$y2),"):  $new";
      }
    }
      }
      else{
	$line{$loc}=$new;
      $score{$loc}=$x-$y;
      }
    }
    else{
push @lines,$_;
    }
  }
  close F;
  open G,(">$file.dedup");
  foreach my $l (@lines){
    if ($l =~ /^\S+:\d+$/){
print G $line{$l};
    }
    else{
print G $l;
    }
  }
  close G;
}
close LS;
