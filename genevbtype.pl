#!/bin/perl
#Created:Jan 21 2017
#Modified:Jan 25 2017
#Q:what atom,bond.angle.. types need to be added? 
#A:it is system-specific, for my system  I add a H-PRO-HSE-ALA-NH2 to the pdb manually(top.0.pdb get). Then I picked up the H on Glu425 to His67 after that I have GluP and Glu , HSD,HSP,HSE in the system(this is done by vmd,top.1.p* get).  Then run ch2lmp to get the data file(top.1.data get). Remove all the things except the types about the H-PRO-HSE-ALA-NH2 added(this is done manually according to a top.4check.data which has no added sequece). At last manually add a hydronium's types(only types).Change some bond parameters into Morse(top.2.data get). The final data file is used as input for this script and for running evb as the initial config
#Q:what kind of atoms (type_name in psf,type_id in lmp.data) invloved given residue names? 
#A:look into the rtf file maybe the best way, so the names provided should be consistent with what appear in the rtf file , or look into the comments in data file is an easy way but may not reliable
#Q:How to know the angles and dih?
#A:read evb.par file one time and then write to evb.type. when reading evb.par .construct an hash of strings, which maps from molecule name to a string "atoms \n bonds\n ...". define global arrays @atoms_in_evbpar,@bonds_in_evbpar..., and def a subrutine to read such a string then construct these global arrays , sub construct_atoms_etc_from_string
#Q:How to change the orders of bond, angle, dihedral,improper in data file so that they are consistent with the evb.par?
#A:generate evb.top by genevbtop.pl first then read it to know what the indexes of  atoms involved in evb(kernel number!=0 and not water, hydronium) and read evb.cfg to now the kernel name. assume the evb.top and evb.par has the same atom order, then we also know what the types of these atoms from evb.par(determined by the third column in evb.top). (data structure: a global hash mapping from evb atom indexes to the kernel number , a hash from indexes to  the types, these two hashes may not be too long because evb molecules in a system is not too many).then loop for all the bonds, when the indexes involved are those evb atoms(check the existency in either the two hashes), look into the  molecule section about that molecule(call construct_atoms_etc_from_string to update the @bonds).loop in @bonds to check which bond is actually involved(do regex twice),get the order of the bond in evb.par(which of the two terminals is larger), if the order in data is the same do nothing, otherwise inverse. loop for angles,dih,impro.
use File::Basename;
use warnings;
$debug=1;$glu_flag=1;$asp_flag=1;
$where_is_this_script=dirname $0;
my %map_from_rtf_atom_type_to_lammps_number_type=qw{};
my %map_from_rtf_bond_type_to_lammps_number_type=qw{};
my %map_from_rtf_angl_type_to_lammps_number_type=qw{};
my %map_from_rtf_dihe_type_to_lammps_number_type=qw{};
my %map_from_rtf_impr_type_to_lammps_number_type=qw{};
my %map_from_lammps_number_type_to_evbpar_bond_type=qw{};
my %map_from_lammps_number_type_to_evbpar_angl_type=qw{};
my %map_from_lammps_number_type_to_evbpar_dihe_type=qw{};
my %map_from_lammps_number_type_to_evbpar_impr_type=qw{};
my @atoms_in_evbpar=qw{},@bonds_in_evbpar=qw{},@angles_in_evbpar=qw{};
my @diherals_in_evbtype=qw{},@impropers_in_evbtype=qw{};#these are just for one kernel
my %map_from_evb_atom_index_to_kernel_name=qw{};
my %map_from_evb_atom_index_to_rtf_type=qw{};
my %kernel_number_to_head_index=qw{};#this works because every evb kernels in the initial config should have different mol_types
sub construct_atoms_etc_from_string {
    my $string=shift;
    unless(defined $string) { die "undefined string\n"; }
    if($string eq "" ) { print "void string in main::construct_atoms_etc_from_string\n"; exit;}
    my @splits=split "\n",$string;
    shift @splits;
    @atoms_in_evbpar=qw{};
    @bonds_in_evbpar=qw{};
    @angles_in_evbpar=qw{};
    @diherals_in_evbtype=qw{};
    @impropers_in_evbtype=qw{};
    #first line tells us the number of atoms,bonds,..
    my @numbers=split ' ',$splits[0];
    @atoms_in_evbpar=splice  @splits,1, $numbers[0];
    @bonds_in_evbpar=splice @splits,1,$numbers[1] ;
    @angles_in_evbpar=splice @splits,1,$numbers[2];
    @dihedrals_in_evbtype=splice @splits,1,$numbers[3];
    @impropers_in_evbtype=splice  @splits ,1,$numbers[4];
}

sub delte_duplicate {
    my $array=shift;
    my %helping_hash=qw{};
    for(@$array) {
        $helping_hash{$_}=1;
    }
    @$array=keys %helping_hash;
}
sub is_equal {
    my $a=shift,$b=shift;
    return 0 unless(keys %$a == keys %$b);
    my %cmp = map { $_ => 1 } keys %$a;
    for my $key (keys %$b) {
        last unless exists $cmp{$key};
        delete $cmp{$key};
    }
    if (%cmp) {
        #print "they don't have the same keys\n";
        return 0;
    } else {
        #print "they have the same keys\n";
        return 1;
    }
}
my $input_count=0;
my @all_the_lines_in_data=qw{};
my @molecules=qw{};
my $mol_global;
my %all_the_atoms_involved=qw{}; #get from rtf based on @molecules\
                             #the hash values are a string containing all the atoms\
                             #involved in the residue
my @all_the_atoms_involved=qw{};
my $map_from_resname_to_kernel_number;
my %map_from_kernel_number_to_molecule_type_section;
print "which data file?\n" if $debug;
while(<STDIN>)
{
    if($input_count==0) {
        chomp;
        if(-e $_)   {$fn_data=$_;}
        else        {$fn_data=$where_is_this_script."/".$_;}
        print "$fn_data\n" if $debug;
        open DATA,"<",$fn_data;
        for(<DATA>) {
            chomp;
            push @all_the_lines_in_data,$_;
        }
        close DATA;
        print "what name of output data?\n" if $debug;
    }
    elsif($input_count==1) {
        print if $debug;
        chomp;
        open DATA,">",$_;
        print "what name of output evb.type?\n" if $debug;
    }
    elsif($input_count==2) {
        print if $debug;
        chomp;
        open TYPE,">",$_;
        print "what name of output bond parameters?\n" if $debug;
    }
    elsif($input_count==3) {
        print if $debug;
        chomp;
        open BOND,">",$_;
        print "what about molecules to declare (name in rtf)?\n" if $debug;
    }
    elsif($input_count==4) {
        chomp;
	$asp_flag=0 unless /ASP/;
	$glu_flag=0 unless /GLU/;
        @molecules=split ' ';
        print "@molecules\n" if $debug;
        print "which rtf file?\n" if $debug;
    }
    elsif($input_count==5) {
        chomp;
        if(-e $_)   {$fn_rtf=$_;}
        else        {$fn_rtf=$where_is_this_script."/".$_;}
        print "$fn_rtf\n" if $debug;
        my $start_flag=0;
        open RTF,"<",$fn_rtf;
        my %atoms_involved_in_this_residue;
        for my $line (<RTF>) {
            chomp;
            my $resname_in_rtf;
            if($line=~m/^RESI/ or $line=~m/^PRES/) {
                $resname_in_rtf=(split ' ' ,$line)[1];
                for $mol (@molecules) {
                    print "$mol\t$resname_in_rtf\n" if $debug>2;
                    if($mol eq $resname_in_rtf) {
                        print "$mol == $resname_in_rtf\n" if $debug>2;
                        $start_flag=1;
                        %atoms_involved_in_this_residue=qw{};
                        $mol_global=$mol;
                        last;
                    }
                }
            }
            if($start_flag==1) {
                if($line=~m/^ATOM/) {
                    print "$line" if $debug>2;
                    $atoms_involved_in_this_residue{(split ' ',$line)[2]}=1;
                }
                if($line=~m/^BOND/) {
                    @atoms_involved_in_this_residue=keys %atoms_involved_in_this_residue; #an array with no duplicates
                    $all_the_atoms_involved{$mol_global}=join " ",@atoms_involved_in_this_residue;
                    $start_flag=0;
                }
            }
        }
################# system specific ############################################
if($glu_flag){
        my $helping_scalar;
        $helping_scalar=$all_the_atoms_involved{"GLU"} if exists $all_the_atoms_involved{"GLU"};                   
        $helping_scalar=$all_the_atoms_involved{"GLU-D"} if exists $all_the_atoms_involved{"GLU-D"};     
        $helping_scalar=$all_the_atoms_involved{"GLU2-D"} if exists $all_the_atoms_involved{"GLU2-D"};
        $helping_scalar=~s/(CC)|(OC)/ /g;                                    
        $all_the_atoms_involved{"GLUP"}=$all_the_atoms_involved{"GLUP"}." ".$helping_scalar if exists  $all_the_atoms_involved{"GLUP"};
        $all_the_atoms_involved{"GLU-P"}=$all_the_atoms_involved{"GLU-P"}." ".$helping_scalar if exists  $all_the_atoms_involved{"GLU-P"};
        $all_the_atoms_involved{"GLU2-P"}=$all_the_atoms_involved{"GLU2-P"}." ".$helping_scalar if exists  $all_the_atoms_involved{"GLU2-P"};
        my @helping_array=qw{};
        @helping_array=split ' ',$all_the_atoms_involved{"GLUP"} if exists  $all_the_atoms_involved{"GLUP"};
        @helping_array=split ' ',$all_the_atoms_involved{"GLU-P"} if exists  $all_the_atoms_involved{"GLU-P"};
        @helping_array=split ' ',$all_the_atoms_involved{"GLU2-P"} if exists  $all_the_atoms_involved{"GLU2-P"};
        my %helping_hash=qw{};
        for(@helping_array) {$helping_hash{$_}=1;}
        @helping_array=keys %helping_hash;
        $all_the_atoms_involved{"GLUP"}=join " ",@helping_array if exists $all_the_atoms_involved{"GLUP"};
        $all_the_atoms_involved{"GLU-P"}=join " ",@helping_array if exists $all_the_atoms_involved{"GLU-P"};
        $all_the_atoms_involved{"GLU2-P"}=join " ",@helping_array if exists $all_the_atoms_involved{"GLU2-P"};
}
if($asp_flag){
        my $helping_scalar=$all_the_atoms_involved{"ASP"};                   
        $helping_scalar=~s/(CC)|(OC)/ /g;                                    
        $all_the_atoms_involved{"ASPP"}=$all_the_atoms_involved{"ASP"}." ".$helping_scalar;
        my @helping_array=split ' ',$all_the_atoms_involved{"ASPP"};
        my %helping_hash=qw{};
        for(@helping_array) {$helping_hash{$_}=1;}
        @helping_array=keys %helping_hash;
        $all_the_atoms_involved{"ASPP"}=join " ",@helping_array;
}
##############################################################################
        for(@molecules) { 
            push @all_the_atoms_involved,(split ' ',$all_the_atoms_involved{$_});
            print "$_:\n$all_the_atoms_involved{$_}\n" if $debug; }
#        print "atom types: @all_the_atoms_involved\n" if $debug;
        close RTF; 
	print "which evb.cfg?\n" if $debug;
    }
    elsif($input_count==6) {
	print if $debug;
	chomp;
	open CFG,"<",$_;
	my $start_declare_molecule=0;
	for(<CFG>) {
		if(/Declare Molecules/) {
			$start_declare_molecule=1;
		}
		if(/Reaction Types/) {
			$start_declare_molecule=0;
		}
		if($start_declare_molecule==1) {
			next if(/^:/ or /^#/ or /^\s*$/);
			my @line_in_cfg=split ' ';
			$map_from_resname_to_kernel_number{$line_in_cfg[1]}=$line_in_cfg[2];
		}
	}
	print "which evb.par?\n" if $debug;
    }
    elsif($input_count==7) {
	print if $debug;
	chomp;
	my $resname_here;
	open PAR,"<",$_;
	my $start_molecule=0,$long_string="";
	for(<PAR>)  {
        chomp;
	    if(/\[molecule_type\.start/) {
		$start_molecule=1;
		/\[.+\..+\.(.+)\]/g;
		$resname_here=$1;
		$long_string="";
	    }
	    if(/\[molecule_type\.end\]/) {
		$start_molecule=0;
	next unless exists $map_from_resname_to_kernel_number{$resname_here};
        $map_from_kernel_number_to_molecule_type_section{$map_from_resname_to_kernel_number{$resname_here}}=$long_string;
        }
	    if($start_molecule==1) {
		next if(/^\s*$/ or /^\[/ or /^:/ or /^\s*#/ or /^\s*:/);
		$long_string=$long_string."\n".$_;	
	   }
	}
	close PAR;
	print %map_from_kernel_number_to_molecule_type_section if $debug>2;
    print "which evb.top?\n" if $debug;
    }
    elsif($input_count==8) {
        print if $debug;
        chomp;
        my $important_flag=0;
        open TOP,"<",$_;
        for(<TOP>) {
            chomp;
            @line_in_top=split ' ';
            unless($line_in_top[1]) {$important_flag=0;next;} 
            my $tmp_resname=$molecules[$line_in_top[1]-1];
            next if($tmp_resname eq "SPCF" or $tmp_resname eq "H3O" or $tmp_resname eq "TIP3" or $tmp_resname eq "H2O"); #we not need to consider the order of bonds.. of water
            if($important_flag==0) {
                $kernel_number_to_head_index{$line_in_top[1]}=$line_in_top[0];
                $important_flag=1;
                print "###@@@ $tmp_resname\n";
            }
           $map_from_evb_atom_index_to_kernel_name{$line_in_top[0]}=$tmp_resname;
        }
        close TOP;
    }
    $input_count++;
    last if($input_count>8) ;
}
###########################################
##           END of STDIN                ## 
###########################################
&delte_duplicate(\@all_the_atoms_involved);
print "\@all_the_atoms_involved:\n@all_the_atoms_involved\n";
my $head_flag=1;$start_masses_flag=0;$start_pair_coeffs_flag=0;
my $start_bond_coeffs_flag=0,$start_atoms_flag=0;
my $start_bonds_flag=0;$start_angle_coeffs_flag=0;$start_angles_flag=0;
my $start_dihedral_coeffs_flag=0;$start_dihedrals_flag=0;
my $start_improper_coeffs_flag=0,$start_improers_flag=0;
for my $line_in_data (@all_the_lines_in_data)
{
    #write to data and bond
    #delete bond coeff of data
    #modify bond coeff , write to bond
    $_=$line_in_data;
    if(/^Masses/) {
        $head_flag=0;$start_masses_flag=1;
    }
    elsif(/^Pair Coeffs/) {
        $start_masses_flag=0;$start_pair_coeffs_flag=1;
    }
    elsif(/^Atoms/) {
        $start_pair_coeffs_flag=0;$start_atoms_flag=1;
    }
    elsif(/^Bond Coeffs/) {
        $start_atoms_flag=0;$start_bond_coeffs_flag=1;
    }
    elsif(/^Bonds/) {
        $start_bond_coeffs_flag=0;$start_bonds_flag=1;
    }
    elsif(/^Angle Coeffs/) {
        $start_bonds_flag=0;$start_angle_coeffs_flag=1;
    }
    elsif(/^Angles/) {
        $start_angle_coeffs_flag=0;$start_angles_flag=1;
    }
    elsif(/^Dihedral Coeffs/) {
        $start_angles_flag=0;$start_dihedral_coeffs_flag=1;
    }
    elsif(/^Dihedrals/) {
        $start_dihedral_coeffs_flag=0;$start_dihedrals_flag=1;
    }
    elsif(/^Improper Coeffs/) {
        $start_dihedrals_flag=0;$start_improper_coeffs_flag=1;
    }
    elsif(/^Impropers/) {
        $start_improper_coeffs_flag=0;$start_improers_flag=1;
    }
    if($head_flag==1) {
        print DATA "$_\n";
    }
    elsif($start_masses_flag==1) {
        print DATA "$_\n";
        for(@all_the_atoms_involved) {
            if($line_in_data=~/#\s+$_$/ or $line_in_data=~/#\s+$_\s+$/) {
                print "find $_ in Masses\n" if $debug;
                my $atom_type=$_;
                if($_ eq "HTH")     {$atom_type="HO";}
                elsif($_ eq "OTH")  {$atom_type="OH";}
                elsif($_ eq "HT")   {$atom_type="HW";}
                elsif($_ eq "OT")   {$atom_type="OW";}
                $map_from_rtf_atom_type_to_lammps_number_type{$atom_type}=(split ' ',$line_in_data)[0];
            }
        }
    } 
    elsif($start_pair_coeffs_flag==1) {
        print DATA "$_\n";
    }
    elsif($start_atoms_flag==1) {
        print DATA "$_\n";
        next if/^\s*$/;
        my $index_in_data=(split ' ')[0];
        my $type=(split ' ')[-1];
        if(exists $map_from_evb_atom_index_to_kernel_name{$index_in_data}) {
            $map_from_evb_atom_index_to_rtf_type{$index_in_data}=$type;
        }
    }
    elsif($start_bond_coeffs_flag==1) {
        #write to par_bond.inp
        next if(/^#/ or /^\s*$/ or /^Bond/);
        my @whole_line=split ' ';
        if(@whole_line<=6) {
            #harmonic
            printf BOND "bond_coeff\t$whole_line[0]\tharmonic %5f %5f\t# $whole_line[4]\t$whole_line[5]\n",$whole_line[1],$whole_line[2];
        }
        else {
            #morse
            print BOND "bond_coeff\t$whole_line[0]\tmorse\t$whole_line[1]\t$whole_line[2]\t$whole_line[3]\t# $whole_line[5]\t$whole_line[6]\n";
        }
        next if @whole_line<5;
        my %bond_name_2compare=qw{};
        my $bond_name_prev="";
        $change_to=$whole_line[-1];
        if($whole_line[-1] eq "HTH")    {$change_to="HO";}
        elsif($whole_line[-1] eq "OTH") {$change_to="OH";}
        $bond_name_2compare{$change_to}=1;
        $change_to=$whole_line[-2];
        if($whole_line[-2] eq "HTH")    {$change_to="HO";}
        elsif($whole_line[-2] eq "OTH") {$change_to="OH";}
        $bond_name_2compare{$change_to}=1;
        for(1..@molecules) {          
            unless(exists $map_from_kernel_number_to_molecule_type_section{$_}) {
               die "$_ not found in hash \$map_from_kernel_number_to_molecule_type_section\n";
            }
            &construct_atoms_etc_from_string
            ($map_from_kernel_number_to_molecule_type_section{$_});
            for(@bonds_in_evbpar) {
                my $bond_name=(split ' ')[-1];
                if($bond_name=~/-/) {@bond_name=split "-",$bond_name;}
                else                {@bond_name=split "_",$bond_name;}
                my %bond_name=qw{};
                if($bond_name[0] eq "BOND") {
                    shift @bond_name;
                }
                for(@bond_name[0..1]) {$bond_name{$_}=1;}
                if(&is_equal(\%bond_name,\%bond_name_2compare)) {
                    if($bond_name ne $bond_name_prev) {
                        print TYPE "#define $bond_name $whole_line[0]\n";
                        $map_from_lammps_number_type_to_evbpar_bond_type
                        {$whole_line[0]}=$bond_name;
                        $bond_name_prev=$bond_name;
                        print "find $bond_name in Bonds\n" if $debug;
                    }
                }
            }
        }   #for(1..@molecules)
    }
    elsif($start_bonds_flag==1) {
        my $modified_line=$_;
        if(/^\s*$/ or /^#/  or /^Bond/) {print DATA "$_\n";next;}
        my @modified_line=split ' ',$modified_line;
        if(exists $map_from_evb_atom_index_to_rtf_type{$modified_line[2]}
            and exists $map_from_evb_atom_index_to_rtf_type{$modified_line[3]}
            and  exists $map_from_lammps_number_type_to_evbpar_bond_type{
             $modified_line[1]}){
            #if all the atoms are evb atoms
            #then it is definitedly an evb bond
            $bond_name_in_evbpar=$map_from_lammps_number_type_to_evbpar_bond_type{
            $modified_line[1]};
            if($map_from_evb_atom_index_to_kernel_name{$modified_line[2]} ne
            $map_from_evb_atom_index_to_kernel_name{$modified_line[3]}) {
                print "not match kernel number, why?\n" if $debug;
            }
            $mol_name=$map_from_evb_atom_index_to_kernel_name{$modified_line[2]};
            unless(exists $map_from_resname_to_kernel_number{$mol_name}) {
               die "$mol_name not found in $map_from_resname_to_kernel_number\n";
            }
            unless(exists $map_from_kernel_number_to_molecule_type_section{$map_from_resname_to_kernel_number{$mol_name}}) {
               die "$map_from_resname_to_kernel_number{$mol_name} not found in $map_from_kernel_number_to_molecule_type_section\n";
            }
                &construct_atoms_etc_from_string(
                $map_from_kernel_number_to_molecule_type_section{
                $map_from_resname_to_kernel_number{$mol_name}});
                for(@bonds_in_evbpar) {
                    if(/$bond_name_in_evbpar/) {
                        @splits=split ' ';
                        $splits[0]+= $kernel_number_to_head_index
                        {$map_from_resname_to_kernel_number{$mol_name}}-1;
                        $splits[1]+= $kernel_number_to_head_index
                        {$map_from_resname_to_kernel_number{$mol_name}}-1;
                        if($splits[0]==$modified_line[3] and $splits[1]==
                        $modified_line[2]) {
                            ($modified_line[2],$modified_line[3])=
                            ($modified_line[3],$modified_line[2]);
                            $modified_line=sprintf "\t%s",join "\t",
                            @modified_line;
                            print "Bond # $modified_line[0] of type $bond_name_in_evbpar is inversed\n" if 
                            $debug;
                        }
                    }
                }
        }
        #else not evb bonds, do nothing
        print DATA "$modified_line\n";
    }
    elsif($start_angle_coeffs_flag==1) {
        print DATA "$_\n";
        my @whole_line=split ' ';
        next if @whole_line<7;
        my $angle_name_prev="";
        my @angle_name_2compare=qw{};
        for(my $iter=-3;$iter<0;$iter++) {
            $change_to=$whole_line[$iter];
            if($whole_line[$iter] eq "HTH")    {$change_to="HO";}
            elsif($whole_line[$iter] eq "OTH") {$change_to="OH";}
            push @angle_name_2compare,$change_to;
        }
        for(1..@molecules) {
           unless(exists $map_from_kernel_number_to_molecule_type_section{$_}) {
              die "$_ not found in $map_from_kernel_number_to_molecule_type_section\n";
           }
            &construct_atoms_etc_from_string
            ($map_from_kernel_number_to_molecule_type_section{$_});
            for(@angles_in_evbpar) {
                my $angle_name=(split ' ',$_)[-1];
               if($angle_name=~/-/) {@angle_name=split "-",$angle_name;}
                else                 {@angle_name=split "_",$angle_name;}
                if($angle_name[0] eq "ANG") {
                    shift @angle_name;
                }
                if(($angle_name[0] eq $angle_name_2compare[0]
                 and $angle_name[2] eq $angle_name_2compare[2] and 
                $angle_name[1] eq $angle_name_2compare[1]) or ($angle_name[0] eq 
                $angle_name_2compare[2] and $angle_name[2] eq 
                $angle_name_2compare[0] and $angle_name[1] eq       
                $angle_name_2compare[1])) 
                {
                    if($angle_name_prev ne $angle_name) {
                        print TYPE "#define $angle_name $whole_line[0]\n";
                        $angle_name_prev=$angle_name;
                        $map_from_lammps_number_type_to_evbpar_angl_type
                        {$whole_line[0]}=$angle_name;
                        print "find $angle_name in Angle Coeffs\n" if $debug;
                    }
                }
            }
        }   #for(1..@molecules)

    }
    elsif($start_angles_flag==1) {
        my $modified_line=$_;
        if(/^\s*$/ or /^#/ or /^Angl/) {print DATA "$_\n";next;} #works
        my @modified_line=split ' ',$modified_line;
#for(@modified_line) {print; print " ";} print "\n";
         if(exists $map_from_evb_atom_index_to_rtf_type{$modified_line[2]}
             and exists $map_from_evb_atom_index_to_rtf_type{$modified_line[3]}
             and exists $map_from_evb_atom_index_to_rtf_type{$modified_line[4]}
	     and exists $map_from_lammps_number_type_to_evbpar_angl_type{
	     $modified_line[1]}){
            #if all the atoms are evb atoms
            #then it is definitedly an evb bond
#            next unless exists $map_from_lammps_number_type_to_evbpar_angl_type{
#            $modified_line[1]};
  #           unless(exists $map_from_lammps_number_type_to_evbpar_angl_type{
 #            $modified_line[1]}) {
#			for(@modified_line) {print; print " ";}
#		}
#		print "\n";
            $angle_name_in_evbpar=$map_from_lammps_number_type_to_evbpar_angl_type{
            $modified_line[1]};
            $mol_name=$map_from_evb_atom_index_to_kernel_name{$modified_line[2]};
            unless(exists $map_from_resname_to_kernel_number{$mol_name}) {
               die "$mol_name not found in $map_from_resname_to_kernel_number\n";
            }
            unless(exists $map_from_kernel_number_to_molecule_type_section{$map_from_resname_to_kernel_number{$mol_name}}) {
               die "$map_from_resname_to_kernel_number{$mol_name} not found in $map_from_kernel_number_to_molecule_type_section\n";
            }
                &construct_atoms_etc_from_string(
                $map_from_kernel_number_to_molecule_type_section{
                $map_from_resname_to_kernel_number{$mol_name}});
                for(@angles_in_evbpar) {
                    if(/$angle_name_in_evbpar/) {
                        @splits=split ' ';
                        $splits[0]+= $kernel_number_to_head_index
                         {$map_from_resname_to_kernel_number{$mol_name}}-1;
                         $splits[2]+= $kernel_number_to_head_index
                         {$map_from_resname_to_kernel_number{$mol_name}}-1;
                        if($splits[0]==$modified_line[4] and $splits[2]==
                        $modified_line[2]) {
#print @splits;print "\n";
#printf "%d\t%d\t%d\n",($splits[0]),($splits[2]),($modified_line[2]-
#$modified_line[4]);
                            ($modified_line[2],$modified_line[4])=
                            ($modified_line[4],$modified_line[2]);
                            $modified_line=sprintf "\t%s",join "\t",
                            @modified_line;
                            print "Angle # $modified_line[0] of type $angle_name_in_evbpar is inversed,WHY?\n"
                            if $debug;
                        }
                    }
                }
        }
        #else not evb bonds, do nothing
        print DATA "$modified_line\n";
    }
    elsif($start_dihedral_coeffs_flag==1) {
        print DATA "$_\n";
        my @whole_line=split ' ';
        next if @whole_line<7;
        my $dihed_name_prev="";
        my @dihed_name_2compare=qw{};
        for(my $iter=-4;$iter<0;$iter++) {
            $change_to=$whole_line[$iter];
            if($whole_line[$iter] eq "HTH")    {$change_to="HO";}
            elsif($whole_line[$iter] eq "OTH") {$change_to="OH";}
            push @dihed_name_2compare,$change_to;
        }
        for(1..@molecules) {
           unless(exists $map_from_kernel_number_to_molecule_type_section{$_}) {
              die "$_ not found in $map_from_kernel_number_to_molecule_type_section\n";
           }
            &construct_atoms_etc_from_string
            ($map_from_kernel_number_to_molecule_type_section{$_});
            for(@dihedrals_in_evbtype) {
                my $dihed_name=(split ' ',$_)[-1];
               if($dihed_name=~/-/) {@dihed_name=split "-",$dihed_name;}
                else                 {@dihed_name=split "_",$dihed_name;}
                if($dihed_name[0] eq "DIH") {
                    shift @dihed_name;
                }
                if(($dihed_name[0] eq $dihed_name_2compare[0]
                 and $dihed_name[1] eq $dihed_name_2compare[1] and 
                $dihed_name[2] eq $dihed_name_2compare[2] and $dihed_name[3] eq 
                $dihed_name_2compare[3]) or ($dihed_name[0] eq 
                $dihed_name_2compare[3] and $dihed_name[1] eq 
                $dihed_name_2compare[2] and $dihed_name[2] eq       
                $dihed_name_2compare[1] and $dihed_name[3] eq 
                $dihed_name_2compare[0])) 
                {
                    if(@dihed_name>4) { #if has suffix 
                        next if($whole_line[2] ne $dihed_name[4]);
                    }
                    if($dihed_name_prev ne $dihed_name) {
                        print TYPE "#define $dihed_name $whole_line[0]\n";
                        $dihed_name_prev=$dihed_name;
                        $map_from_lammps_number_type_to_evbpar_dihe_type
                        {$whole_line[0]}=$dihed_name;
                        print "find $dihed_name in Dihedral Coeffs\n" if $debug;
                    }
                }
            }
        }   #for(1..@molecules)

    }
    elsif($start_dihedrals_flag==1) {
       my $modified_line=$_;
       if(/^\s*$/ or /^#/ or /^Dihed/) {print DATA "$_\n";next;}
        my @modified_line=split ' ',$modified_line;
        if(exists $map_from_evb_atom_index_to_rtf_type{$modified_line[2]}
            and exists $map_from_evb_atom_index_to_rtf_type{$modified_line[3]}
            and exists $map_from_evb_atom_index_to_rtf_type{$modified_line[4]}  
            and exists $map_from_evb_atom_index_to_rtf_type{$modified_line[5]}
	    and exists $map_from_lammps_number_type_to_evbpar_dihe_type{ 
            $modified_line[1]}){
            #if all the atoms are evb atoms
            #then it is definitedly an evb bond
            $dihedral_name_in_evbpar=
            $map_from_lammps_number_type_to_evbpar_dihe_type{
            $modified_line[1]};
            $mol_name=$map_from_evb_atom_index_to_kernel_name{$modified_line[2]};
            unless(exists $map_from_resname_to_kernel_number{$mol_name}) {
               die "$mol_name not found in $map_from_resname_to_kernel_number\n";
            }
            unless(exists $map_from_kernel_number_to_molecule_type_section{$map_from_resname_to_kernel_number{$mol_name}}) {
               die "$map_from_resname_to_kernel_number{$mol_name} not found in $map_from_kernel_number_to_molecule_type_section\n";
            }
                &construct_atoms_etc_from_string(
                $map_from_kernel_number_to_molecule_type_section{
                $map_from_resname_to_kernel_number{$mol_name}});
#print "\@dih= @dihedrals_in_evbtype\n"; 
#printf "%s\n",$map_from_resname_to_kernel_number{$mol_name};
                for(@dihedrals_in_evbtype) {
                    if(/$dihedral_name_in_evbpar/) {
                        @splits=split ' ';
#print $kernel_number_to_head_index
                        #{$map_from_resname_to_kernel_number{$mol_name}}."\n";
                        $splits[0]+= $kernel_number_to_head_index
                          {$map_from_resname_to_kernel_number{$mol_name}}-1;
                          $splits[3]+= $kernel_number_to_head_index
                          {$map_from_resname_to_kernel_number{$mol_name}}-1;
#print "$splits[0]  $splits[3]\n";
                         if($splits[0]==$modified_line[5] and 
                            $splits[3]==$modified_line[2] ) {
                            ($modified_line[2],$modified_line[5])=
                            ($modified_line[5],$modified_line[2]);
                            ($modified_line[3],$modified_line[4])=
                            ($modified_line[4],$modified_line[3]);
                            $modified_line=sprintf "\t%s",join "\t",
                            @modified_line;
                            print "Dihedral # $modified_line[0] of type $dihedral_name_in_evbpar is inversed\n"
                            if $debug;
                        }
                    }
                }
        }
        #else not evb bonds, do nothing
        print DATA "$modified_line\n";
    }
    elsif($start_improper_coeffs_flag==1) {
        print DATA "$_\n";
        my @whole_line=split ' ';
        next if @whole_line<7;
        my $impro_name_prev="";
        my @impro_name_2compare=qw{};
        for(my $iter=-4;$iter<0;$iter++) {
            $change_to=$whole_line[$iter];
            if($whole_line[$iter] eq "HTH")    {$change_to="HO";}
            elsif($whole_line[$iter] eq "OTH") {$change_to="OH";}
            push @impro_name_2compare,$change_to;
        }
#for(@impro_name_2compare) {print ; print " ";} print "\n";
        for(1..@molecules) {
           unless(exists $map_from_kernel_number_to_molecule_type_section{$_}) {
              die "$_ not found in $map_from_kernel_number_to_molecule_type_section\n";
           }
            &construct_atoms_etc_from_string
            ($map_from_kernel_number_to_molecule_type_section{$_});
            for(@impropers_in_evbtype) {
                my $impro_name=(split ' ',$_)[-1];
               if($impro_name=~/-/) {@impro_name=split "-",$impro_name;}
                else                 {@impro_name=split "_",$impro_name;}
                if($impro_name[0] eq "IMP") {
                    shift @impro_name;
                }
                if(($impro_name[0] eq $impro_name_2compare[0]
                 and $impro_name[1] eq $impro_name_2compare[1] and 
                $impro_name[2] eq $impro_name_2compare[2] and $impro_name[3] eq 
                $impro_name_2compare[3]) or ($impro_name[0] eq 
                $impro_name_2compare[3] and $impro_name[1] eq 
                $impro_name_2compare[2] and $impro_name[2] eq       
                $impro_name_2compare[1] and $impro_name[3] eq 
                $impro_name_2compare[0])) 
                {
                    if(@impro_name>4) { #if has suffix 
                        next if($whole_line[2] ne $impro_name[4]);
                    }
                    if($impro_name_prev ne $impro_name) {
                        print TYPE "#define $impro_name $whole_line[0]\n";
                        $impro_name_prev=$impro_name;
                        $map_from_lammps_number_type_to_evbpar_impr_type
                        {$whole_line[0]}=$impro_name;
                        print "find $impro_name in Improral Coeffs\n" if $debug;
                    }
                }
            }
        }   #for(1..@molecules)

    }
    elsif($start_improers_flag==1) {
       my $modified_line=$_;
       if(/^\s*$/ or /^#/ or /^Impro/) {print DATA "$_\n";next;}
        my @modified_line=split ' ',$modified_line;
        if(exists $map_from_evb_atom_index_to_rtf_type{$modified_line[2]}
            and exists $map_from_evb_atom_index_to_rtf_type{$modified_line[3]}
            and exists $map_from_evb_atom_index_to_rtf_type{$modified_line[4]}  
            and exists $map_from_evb_atom_index_to_rtf_type{$modified_line[5]}
	    and exists $map_from_lammps_number_type_to_evbpar_impr_type{
            $modified_line[1]}){
            #if all the atoms are evb atoms
            #then it is definitedly an evb bond
            $improper_name_in_evbpar=
            $map_from_lammps_number_type_to_evbpar_impr_type{
            $modified_line[1]};
            $mol_name=$map_from_evb_atom_index_to_kernel_name{$modified_line[2]};
            unless(exists $map_from_resname_to_kernel_number{$mol_name}) {
               die "$mol_name not found in $map_from_resname_to_kernel_number\n";
            }
            unless(exists $map_from_kernel_number_to_molecule_type_section{$map_from_resname_to_kernel_number{$mol_name}}) {
               die "$map_from_resname_to_kernel_number{$mol_name} not found in $map_from_kernel_number_to_molecule_type_section\n";
            }
                &construct_atoms_etc_from_string(
                $map_from_kernel_number_to_molecule_type_section{
                $map_from_resname_to_kernel_number{$mol_name}});
                for(@impropers_in_evbtype) {
                    if(/$improper_name_in_evbpar/) {
                        @splits=split ' ';
                        $splits[0]+= $kernel_number_to_head_index
                          {$map_from_resname_to_kernel_number{$mol_name}}-1;
                          $splits[3]+= $kernel_number_to_head_index
                          {$map_from_resname_to_kernel_number{$mol_name}}-1;
                         if($splits[0]==$modified_line[5] and 
                            $splits[3]==$modified_line[2] ) {
                            ($modified_line[2],$modified_line[5])=
                            ($modified_line[5],$modified_line[2]);
                            ($modified_line[3],$modified_line[4])=
                            ($modified_line[4],$modified_line[3]);
                            $modified_line=sprintf "\t%s",join "\t",
                            @modified_line;
                            print "Improper # $modified_line[0] of type $improper_name_in_evbpar is inversed\n"
                            if $debug;
                        }
                    }
                }
        }
        #else not evb bonds, do nothing
        print DATA "$modified_line\n";
    }
}
close DATA;
for(keys %map_from_rtf_atom_type_to_lammps_number_type) {
    print TYPE "#define $_ $map_from_rtf_atom_type_to_lammps_number_type{$_}\n";
}
close BOND;
###############################system specific######################
#print TYPE "#define HNd 67\n";
####################################################################
close TYPE;
#$dooooo=5;
#&construct_atoms_etc_from_string($map_from_kernel_number_to_molecule_type_section{$dooooo});
#for(@impropers_in_evbtype) {print ; print " ";}
#print %map_from_lammps_number_type_to_evbpar_bond_type;
#print %map_from_evb_atom_index_to_kernel_name;
#print "\n";
#print %map_from_evb_atom_index_to_rtf_type;
#print "\n";
#print "\n@atoms_in_evbpar\n\n@bonds_in_evbpar\n\n@angles_in_evbpar\n\n@diherals_in_evbtype\n\n@impropers_in_evbtype\n";
#print %map_from_lammps_number_type_to_evbpar_bond_type;
#print "\n";
#print %map_from_lammps_number_type_to_evbpar_angl_type;
#print "\n";
#print %map_from_lammps_number_type_to_evbpar_dihe_type;
#print "\n";
#print %map_from_lammps_number_type_to_evbpar_impr_type;
#print "\n";
#printf "%s", join " " ,keys %map_from_evb_atom_index_to_rtf_type;print "\n";
#for(keys %kernel_number_to_head_index) {print "$_  $kernel_number_to_head_index{$_} \n";}
#for(keys %map_from_evb_atom_index_to_kernel_name) {print "$_  $map_from_evb_atom_index_to_kernel_name{$_} \n";}
