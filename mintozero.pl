use strict;
my @input_files = @ARGV;
my $min = 999999;
while(<>) {
   next if /#/;
   my @F = split ' ';
   $min = $F[1] if $F[1] < $min;
}
@ARGV = @input_files;
while(<>) {
   next if /#/;
   my @F = split ' ';
   printf "$F[0] %f\n", $F[1] - $min;
}