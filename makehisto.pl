# max = min + width*max_nbin
# max and min must be suitable for the real ranges !!
$width=0.01;
$col=19;
$min=1.85;
$max_nbin=145;
@histo=qw{};
push @histo,0.0 for(0..$max_nbin);
while(<>) {
    next if /^#/;
    @F=split ' ';
    $histo[int(($F[$col]-$min)/$width)]+=1;
}
for(0..$max_nbin) {
    printf "%lf %d\n",$_*$width+$min+$width/2.0,$histo[$_];
}
