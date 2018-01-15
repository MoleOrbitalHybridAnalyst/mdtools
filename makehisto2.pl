# perl makehisto2.pl nbins column_index file1 file2 ...
# more intelligent version of makehisto.pl
use strict;
#### hard coded ####
my $ratio = 0.2; # extend the histo range by 0.2*actual range
####################
my $max_nbin = shift @ARGV;
my $col = shift @ARGV;
my @input_files = @ARGV;
# first loop to find ranges
my $count = 0; my $min; my $max;
while(<>) {
   next if /^#/ or /^\s*$/;
   my @F=split ' ';
   my $value = $F[$col];
   if($count==0) {
      $max = $min = $value;
      $count++; next;
   }
   $max = $value if $value>$max;
   $min = $value if $value<$min;
}
my $extend = ($max - $min) * $ratio/2.0;
$min -= $extend; $max += $extend;
my $width = ($max - $min)/$max_nbin;
#@@@ print "$max $min $max_nbin\n";
# second loop to make histo
my @histo=qw{};
push @histo,0.0 for(0..$max_nbin);
@ARGV = @input_files; #rewind input files
while(<>) {
   next if /^#/ or /^\s*$/;
   my @F=split ' ';
   $histo[int(($F[$col]-$min)/$width)]+=1;
}
for(0..$max_nbin) {
    printf "%lf %d\n",$_*$width+$min+$width/2.0,$histo[$_];
}
