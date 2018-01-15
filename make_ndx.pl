# perl make_ndx file group_name1 group_name2 ...
use strict;

my $num_groups = $#ARGV;
my $input_file = shift @ARGV;
my @group_names = @ARGV;
@ARGV=qw{}; push @ARGV, $input_file;
my $line_count = 0;
while(<>) {
   next if /[#@!]/;
   $line_count++;
   last if($line_count > $num_groups);
   chomp;
   my @F = split ' ';
   printf "[ %s ]\n", $group_names[$line_count - 1];
   my $ele_count = 0;
   for(@F) {
      my $c = " ";
      $c = "\n" if ($ele_count + 1) % 8 == 0;
      printf "%d%s",$_,$c;
      $ele_count++;
   }
   print "\n\n";
}
