#!/usr/bin/perl
# perl make_pandas_happy.pl plumed_colvar 
# convert plumed colvar to csv
use strict;
my $header_count = 0;
while(<>) {
   chomp;
   if(/^#/) {
      next if $header_count != 0;
      $header_count++;
      my @F = split ' ';
      for my $f (@F[2..$#F-1]) {
         print "$f,";
      }
      print "$F[$#F]\n"; next;
   }
   my @F = split ' ';
   for my $f (@F[0..$#F-1]) {
      print "$f,";
   }
   print "$F[$#F]\n";
}

