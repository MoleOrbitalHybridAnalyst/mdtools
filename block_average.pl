# perl block_average.pl nblocks histo_min histo_max nbins 
# should have ../metadatafile.full
use strict;
my $nblocks = $ARGV[0];
my @targets = qw{};
my @kappas = qw{};
my @dir_names = qw{};
open IN, "<", "../metadatafile.full";
while(<IN>) {
   chomp;
   my @F=split ' ';
   $_=$F[0];
   push @targets, $F[1];
   push @kappas, $F[2];
   /\/([m0-9\-\._]+)\/(COLVAR|tail_colvar)/;
   my $dir_name = $1;
   push @dir_names, $dir_name;
   system "mkdir $dir_name";
   chdir $dir_name;
   my $nline = int(`wc -l ../../$F[0]|cut -d' ' -f1` / $nblocks);
   system "split -l$nline ../../$F[0] -da2";
   chdir "../";
}
my $histo_min=$ARGV[1];
my $histo_max=$ARGV[2];
my $nbins=$ARGV[3];
close IN;
for my $id (0..($nblocks-1)) {
   open OUT, ">", "metadatafile".$id;
   my $i = 0;
   for(@dir_names) {
      print OUT "$_/x";
      printf OUT "%.2d\t", $id;
      printf OUT "%f\t", $targets[$i];
      printf OUT "%f\n", $kappas[$i];
      $i ++;
   }
   close OUT;
   system "wham $histo_min $histo_max $nbins 1e-11 310 0 metadatafile${id} pmf${id}.dat > wham${id}.log";
}
