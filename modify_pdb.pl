#editconf gives wrong atom type as shown in step8_3.pdb
#need to fix the wrong atom type
#and need to give correct segname
$num1=7953;
$num2=23328;
$num3=28535;
$num4=56477;
$num5=61499;
$num6=61536;
while(<>)
{
	chomp;
	if(/^ATOM/)
	{
		@s1=split ' ';
		$ss=$s1[2]; #atom type
		if($ss=~ m/^\d/){
			@s2=split '',$ss;
			my $sw=$s2[0];
			for($i=1;$i<@s2;$i++)
			{
				$s2[$i-1]=$s2[$i];
			}
			$s2[-1]=$sw;
			my $s2="";
			for(@s2)
			{
				$s2=join("",$s2,$_);
			}
			s/$ss/$s2/e;
		}
		s/\s*$//;
		$_=$_."      ";
		if($s1[1]<$num1) { $_=$_."PROA"; }
		elsif($s1[1]<$num2) { $_=$_."POPE"; }
        elsif($s1[1]<$num3) { $_=$_."POPG"; }
		elsif($s1[1]<$num4) { $_=$_."W1";}
		elsif($s1[1]<$num5) { $_=$_."W2";}
        elsif($s1[1]<$num6) { $_=$_."SOD";}
		else { $_=$_."CLA"; }
	}
	print $_."\n";
}
