while(<>) {
   chomp;
   next if /^#/;
   my @F = split ' ';
   if($F[1]=="inf") {
      $F[1] = 50;
   }
   print "$F[0] $F[1]\n";
}
