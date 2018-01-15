# perl force_diff.pl col_index1 col_fx1 col_fy1 col_fz1 col_index2 col_fx2 col_fy2 col_fz2 force_file1 force_file2
use strict;
my $col_index1 = shift @ARGV;
my $col_fx1 = shift @ARGV;
my $col_fy1 = shift @ARGV;
my $col_fz1 = shift @ARGV;
my $col_index2 = shift @ARGV;
my $col_fx2 = shift @ARGV;
my $col_fy2 = shift @ARGV;
my $col_fz2 = shift @ARGV;
my $force_file1 = shift @ARGV;
my $force_file2 = shift @ARGV;

my @indexes = qw{};
my @fx1 = qw{};
my @fy1 = qw{};
my @fz1 = qw{};

open F1, "<", $force_file1;
while(<F1>) {
   next if /^#!@$/;
   chomp;
   my @F = split ' ';
   push @indexes, $F[$col_index1];
   push @fx1, $F[$col_fx1];
   push @fy1, $F[$col_fy1];
   push @fz1, $F[$col_fz1];
}
close F1;

open F2, "<", $force_file2;
while(<F2>) {
   next if /^#!@$/;
   chomp;
   my @F = split ' ';
   my $i = 0;
   for(@indexes) {
      if($_ == $F[$col_index2]) {
         printf "%d %f %f %f\n", 
            $_, 
            $fx1[$i]-$F[$col_fx2],
            $fy1[$i]-$F[$col_fy2],
            $fz1[$i]-$F[$col_fz2];
         last;
      }
      $i++;
   }
}
close F2;
