use strict;
use POSIX qw{floor};

my $prdx = 19.388697220179999;
my $prdy = 19.388697220179999;
my $prdz = 19.388697220179999;

sub pbc() {
   my $x = shift;
   my $refx = shift;
   my $b = shift;
   my $k = floor(($x - $refx)/$b + 0.5);
   return ($x - $refx) - $k * $b + $refx;
}

my $stat = 0;
my $step = 0;
my $pdx; my $pdy; my $pdz;
my $dx; my $dy; my $dz;
$pdx = $pdy = $pdz = 0.0;
$dx = $dy = $dz = 0.0;
while(<>) {
   if(/DIPOLE/) {
      $stat ++;
      next;
   }
   if($stat == 1) {
      $pdx = $dx;
      $pdy = $dy;
      $pdz = $dz;
      my @F = split ' ';
      $dx = $F[0]; $dy = $F[1]; $dz = $F[2];
      if($step) {
         $dx = &pbc($dx, $pdx, $prdx);
         $dy = &pbc($dy, $pdy, $prdy);
         $dz = &pbc($dz, $pdz, $prdz);
      }
      print "$step $dx $dy $dz\n";
      $stat = 0;
      $step ++;
   }
}
