#!/usr/bin/perl -w
use strict;
#N.B. this requires 23andme_similarity_matrix.pl to be run first (SLOW)
#look at max_similarity_matrix.tab to choose reasonable cutoff for max similarity and minimum similarity
#min similraity files are removed, but max similarity files you keep one of them
#check for and remove checkpoint files (tempfiles) if running first time
my $percentid;
my (@files,@size);
my %exclude;
my $simfile='similarity_matrix.tab';
my $maxsimfile='max_similarity_matrix.tab';
my $usage="23andme_similarity_filter.pl <max percentid> <#min percentid>";
die "Usage: $usage\n" unless (@ARGV >= 1);
my $filter = $ARGV[0];
my $minfilter = 0; if (defined($ARGV[1])){$minfilter = $ARGV[1];}
my %list;

 open L,("filelist.list");
 while (<L>){
   if (/^(\d+)\t(\S.*)$/){
     my $file=$2;chomp($file);push @files,$1;
     $list{$1}=$file;
   }
 }
 close L;

open T,($simfile);
while (<T>){
  if (/^(\d+)\t(\d+)\t(\S+)/){
    unless (exists($exclude{$1}) or exists($exclude{$2})){
    if (100*$3 > $filter){
      $exclude{$2}=1;
      print STDERR "$2 excluded because of ",abs(100*$3),"% similarity to $1\n";
  }
  }
  }
}
close T;

open T,($maxsimfile);
while (<T>){
  if (/^(\S+)\t(\d+)$/){
    if (100*$1 <= $minfilter){
$exclude{$2}=1;
      print STDERR "$2 excluded because of ",abs(100*$1),"% max similarity to something else\n";
    }
  }
}
close T;
   

foreach my $file (@files){
  unless (exists($exclude{$file})){
print "$file\t$list{$file}\n";
  }
}
