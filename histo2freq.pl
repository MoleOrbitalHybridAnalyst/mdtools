use strict;
#my @input_files = @ARGV;
my $total = 0;
my @x = qw{};
my @y = qw{};
while(<>) {
   next if /#/;
   my @F = split ' ';
   push @y, $F[1];
   push @x, $F[0];
}
#@ARGV = @input_files;
#while(<>) {
#   next if /#/;
#   my @F = split ' ';
#   printf "$F[0] %f\n", $F[1] / $total;
#}
for(@y) {
   $total += $_;
}
for(0..$#y) {
   printf "$x[$_] %f\n", $y[$_] / $total;
}
