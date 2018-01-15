while(<>) {
   chomp;
   next if /^#/;
   my @F = split ' ';
   print "$F[0] $F[1]\n";
}
