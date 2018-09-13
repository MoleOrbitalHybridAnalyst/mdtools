use strict;
my @input_files = @ARGV;
my $total = 0;
while(<>) {
   next if /#/;
   my @F = split ' ';
   $total += $F[1];
}
@ARGV = @input_files;
while(<>) {
   next if /#/;
   my @F = split ' ';
   printf "$F[0] %f\n", $F[1] / $total;
}
