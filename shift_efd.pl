use strict;
my $offset = 0;
#my $should_shift = 0;
my $efd_value = 0;
while(<>) {
   if(/#/) {
#      $should_shift = 1;
      $offset = $efd_value;
      next;
   }
   chomp;
   my @F = split ' ';
   $efd_value = $offset + $F[1];
   print "$F[0] $efd_value\n";
#   if($should_shift) 
}
