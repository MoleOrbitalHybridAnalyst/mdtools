# pbc(x,box)
sub pbc() {
   my $x = shift;
   my $box = shift;
   if($x>$box/2.0) {
      $x -= $box;
   }
   elsif($x<-$box/2.0) {
      $x += $box;
   }
   else {
      return $x;
   }
   $x = &pbc($x,$box);
}
1;
