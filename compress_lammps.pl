my $offset=9; # # of lines = $natom + $offset
my $natom=$ARGV[0];
my $ratio=$ARGV[1];
my $nline=$natom + $offset;
my $count = 0;
my $subcount = 0;
my $isprint = 0;
@ARGV=@ARGV[2..@ARGV-1];
while(<>) {
    if($count % ($ratio * $nline) == 0) {
        $isprint = 1;
        $subcount = 0;
    }
    if($isprint) {
        print;
        $subcount++;
    }
    $isprint = 0 if $subcount == $nline; 
    $count++;
}
