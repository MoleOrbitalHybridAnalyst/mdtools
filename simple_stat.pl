#!/bin/perl

use strict;
use Getopt::Long;

my $cols = '0';
GetOptions ("col=s" => \$cols);
my @cols = split /,/, $cols;

my $i = 0;
my @values = qw{};
while(<>) {
   next if /#/;
   chomp;
   $i++;
   my @F = split ' ';
   for my $icol (0..$#cols) {
      push @{ $values[$i] }, $F[$cols[$icol]];
   }
}

for my $icol (0..$#cols) {
   my $max = $values[1][$icol];
   my $min = $values[1][$icol];
   my $argmax = 1;
   my $argmin = 1;
   my $absmax = abs($values[1][$icol]);
   my $absmin = abs($values[1][$icol]);
   my $argabsmax = 1;
   my $argabsmin = 1;
   my $ave = 0;
   for my $j (1..$i) {
      my $value = $values[$j][$icol];
      $ave += $value;
      if($value > $max) {
         $max = $value; $argmax = $j;
      }
      if($value < $min) {
         $min = $value; $argmin = $j;
      }
      if(abs($value) > $absmax) {
         $absmax = abs($value); $argabsmax = $j;
      }
      if(abs($value) < $absmin) {
         $absmin = abs($value); $argabsmin = $j;
      }
   }
   $ave /= $i;
   my $var = 0.0;
   for my $j (1..$i) {
      my $value = $values[$j][$icol];
      $var += ($value - $ave) ** 2;
   }
   $var /= $i;
   print "column $icol:\n";
   print "max = $max at line $argmax\n";
   print "min = $min at line $argmin\n";
   print "absmax = $absmax at line $argabsmax\n";
   print "absmin = $absmin at line $argabsmin\n";
   print "average = $ave\n";
   print "variance (dof = 0)  = $var\n";
}
