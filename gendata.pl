#perl gendata.pl check.data xxx.1.data xxx.2.data [resid1 resid2 ...]
#modify xxx.1.data to xxx.2.data
    #remove extra atoms bonds angles dihedrals impropers according to check.data
    #change # of atoms, # of bonds, ... accordingly
    #add HNd to the end of atom type list
        #masses; pair coeffs
        #change atom type of H in residue resid1, resid2, ...  to the type of HNd; do not need to change bonds, angles, ... because thoses have their own types
    #modify bond coeffs of the reactive residues to use morse potential
        #H NR3 , HTH OTH

use strict;
use warnings;

my $natom=0;my $nbond=0;my $nangle=0;my $ndihedral=0;my $nimproper=0; 
my $ncrossterms=0;
open CHECK, "<", $ARGV[0];
open DATA2, ">", $ARGV[2];
my $time=localtime;
print DATA2 "Created by gendata.pl on $time\n\n";
for(<CHECK>) {
    chomp;
    if(/^\s+[0-9]+\s+atoms/) {
        $natom=(split ' ')[0];
        print DATA2 $_."\n";
    }
    elsif(/^\s+[0-9]+\s+bonds/) {
        $nbond=(split ' ')[0];
        print DATA2 $_."\n";
    }
    elsif(/^\s+[0-9]+\s+angles/) {
        $nangle=(split ' ')[0];
        print DATA2 $_."\n";
    }
    elsif(/^\s+[0-9]+\s+dihedrals/) {
        $ndihedral=(split ' ')[0];
        print DATA2 $_."\n";
    }
    elsif(/^\s+[0-9]+\s+impropers/) {
        $nimproper=(split ' ')[0];
        print DATA2 $_."\n";
    }
    elsif(/^\s+[0-9]+\s+crossterms/) {
        $ncrossterms=(split ' ')[0];
        print DATA2 $_."\n";
        last;
    }
}
close CHECK;
print DATA2 "\n";

open DATA1, "<", $ARGV[1];
my $headflag=0;my $Masses_flag=0; my $Pair_Coeffs_flag=0; my $Atoms_flag=0; my $Bond_Coeffs_flag=0; my $Bonds_flag=0; my $Angle_Coeffs_flag=0; my $Angles_flag=0; my $Dihedral_Coeffs_flag=0; my $Dihedrals_flag=0; my $Improper_Coeffs_flag=0; my $Impropers_flag=0; my $CMAP_flag=0;
my $natomtypes=0;
my $pair_coeffs_for_h="";
my $type_for_h=0;
my $count=-1;
my $HNd_count=0;
my $flag_mod=0;
for(<DATA1>) {
    chomp;
    if(/^\s+[0-9]+\s+atom types/) {
        $headflag=1;
    }
    elsif(/^Masses/) {
        $count=-1;
        $Masses_flag=1;$headflag=0;
    }
    elsif(/^Pair Coeffs/) {
        $count=-1;
        $Pair_Coeffs_flag=1;$Masses_flag=0;
    }
    elsif(/^Atoms/) {
        $count=-1;
        $Atoms_flag=1;$Pair_Coeffs_flag=0;
    }
    elsif(/^Bond Coeffs/) {
        $count=-1;
        $Bond_Coeffs_flag=1; $Atoms_flag=0;
    }
    elsif(/^Bonds/) {
        $count=-1;
        $Bonds_flag=1; $Bond_Coeffs_flag=0;
    }
    elsif(/^Angle Coeffs/) {
        $count=-1;
        $Angle_Coeffs_flag=1; $Bonds_flag=0;
    }
    elsif(/^Angles/) {
        $count=-1;
        $Angles_flag=1; $Angle_Coeffs_flag=0;
    }
    elsif(/^Dihedral Coeffs/) {
        $count=-1;
        $Dihedral_Coeffs_flag=1; $Angles_flag=0;
    }
    elsif(/^Dihedrals/) {
        $count=-1;
        $Dihedrals_flag=1; $Dihedral_Coeffs_flag=0;
    }
    elsif(/^Improper Coeffs/) {
        $count=-1;
        $Improper_Coeffs_flag=1; $Dihedrals_flag=0;
    }
    elsif(/^Impropers/) {
        $count=-1;
        $Impropers_flag=1; $Improper_Coeffs_flag=0;
    }
    elsif(/^CMAP/) {
        $count=-1;
        $CMAP_flag=1; $Impropers_flag=0;
    }
    if($headflag==1) {
        if(/^\s+[0-9]+\s+atom types/) {
            if(@ARGV>3) {
                $natomtypes=(split ' ')[0]+1;
            }
            else {
                $natomtypes=(split ' ')[0];
            }
            printf DATA2 "%12d  atom types\n",$natomtypes;
        }
        else {
            print DATA2 $_."\n";
        }
    }
    elsif($Masses_flag==1) {
        print DATA2 $_."\n";
        $count++ unless /^\s*$/;
        if($count==$natomtypes-1 and @ARGV>3) {
            printf DATA2 "%8d %10.7g%s\n",$natomtypes,1.008,"  # HNd";
            $count++;
        }
    }
    elsif($Pair_Coeffs_flag==1) {
        print DATA2 $_."\n";
        $count++ unless /^\s*$/;
        if(/#\sH\s*$/) {
             $pair_coeffs_for_h=$_;
        }
        if($count==$natomtypes-1 and @ARGV>3) {
            my @tmp=split ' ',$pair_coeffs_for_h;
            $type_for_h=$tmp[0];
            printf DATA2 "%8d %10.7g %10.7g %10.7g %10.7g%s\n",$natomtypes,$tmp[1],$tmp[2],$tmp[3],$tmp[4],"  # HNd";
            $count++;
        }
    }
    elsif($Atoms_flag==1) {
        $flag_mod=0;
        $count++ unless /^\s*$/;
        next if $count>$natom;
        my @tmp=split ' ';
        if(@tmp>1 and @ARGV>3) {
            foreach my $resid (@ARGV[3..@ARGV-1]) {
                if($resid==$tmp[1]) {
                    if($type_for_h==$tmp[2]) {
                        $HNd_count++;
                        if($HNd_count>1) {
                            printf DATA2 "%8d %7d %5d %9.6g %11.8g %11.8g %11.8g%s\n",$tmp[0],$tmp[1],$natomtypes,$tmp[3],$tmp[4],$tmp[5],$tmp[6],"  # H";
                            print "atom # $tmp[0] is changed into type $natomtypes\n";
                            $HNd_count=0;
                            $flag_mod=1;
                        }
                    }
                }
            }
        }
        print DATA2 $_."\n" if $flag_mod==0;
        print DATA2 "\n" if $count==$natom;
    }
    elsif($Bond_Coeffs_flag==1) {
        $flag_mod=0;
        $count++ unless /^\s*$/;
        if(/H\s+NR3/) {
            print DATA2 "$count\t135.07\t2.06\t1.01\t# H    NR3\n";
            $flag_mod=1;
        }
        if(/HTH\s+OTH/) {
            print DATA2 "$count\t79.0864\t2.0834\t0.98\t# HTH  OTH\n";
            $flag_mod=1;
        }
        if(/H\s+OH1/) {
            print DATA2 "$count\t143.003\t1.80\t0.975\t# H   OH1\n";
            $flag_mod=1;
        }
        print DATA2 $_."\n" if $flag_mod==0;
    }
    elsif($Bonds_flag==1) {
        $count++ unless /^\s*$/;
        next if $count>$nbond;
        print DATA2 $_."\n";
        if($count==$nbond) {
           print DATA2 "\n";
           $count++;
           next;
        }
    }
    elsif($Angle_Coeffs_flag==1) {
        $count++ unless /^\s*$/;
        print DATA2 $_."\n";
    }
    elsif($Angles_flag==1) {
        $count++ unless /^\s*$/;
        next if $count>$nangle;
        print DATA2 $_."\n";
        if($count==$nangle) {
           print DATA2 "\n";
           $count++;
           next;
        }
    }
    elsif($Dihedral_Coeffs_flag==1) {
        $count++ unless /^\s*$/;
        print DATA2 $_."\n";
    }
    elsif($Dihedrals_flag==1) {
        $count++ unless /^\s*$/;
        next if $count>$ndihedral;
        print DATA2 $_."\n";
        if($count==$ndihedral) {
           print DATA2 "\n";
           $count++;
           next;
        }
    }
    elsif($Improper_Coeffs_flag==1) {
        $count++ unless /^\s*$/;
        print DATA2 $_."\n";
    }
    elsif($Impropers_flag==1) {
        $count++ unless /^\s*$/;
        next if $count>$nimproper;
        print DATA2 $_."\n";
        if($count==$nimproper) {
           print DATA2 "\n";
           $count++;
           next;
        }
    }
    elsif($CMAP_flag==1) {
        $count++ unless /^\s*$/;
        next if $count>$ncrossterms;
        print DATA2 $_."\n";
        if($count==$ncrossterms) {
           print DATA2 "\n";
           $count++;
           next;
        }
    }
}
close DATA1;
close DATA2;
