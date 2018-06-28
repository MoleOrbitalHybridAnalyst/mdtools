#!/bin/perl
#read top.0.pdb to get every atom's resname resid segname, ask for inputs of totoal atom numbers in the system (61544 in this case), all the resname resid segname to be treated evb-ly and then check for every atom in the pdb . also read evb.cfg to get the kernel number of the residue.
#####MAYBE BUGGY IF TWO EVB RESIDUES APPEAR ADJECENTLY IN PDB FILE
use warnings;
$debug=1;
$input_count=0;
my $total_number_atoms;
my %map_from_resname_to_kernerl_number;
my @resname_resid_segname=qw{};
my @resnames=qw{},@resids=qw{},@segnames=qw{};
print "what is the totoal number of atoms in the system?\n" if $debug;
while(<STDIN>) {
    if($input_count==0) {
        print if $debug;
        chomp;
        $total_number_atoms=$_;
        print "which evb.cfg?\n";
    }
    elsif($input_count==1) {
        print if $debug;
        chomp;
        open CFG,"<",$_;
        my $start_declare_molecule=0;
        for(<CFG>) {
            if(/Declare Molecules/) {
                $start_declare_molecule=1;
            }
            if($start_declare_molecule==1) {
                $start_declare_molecule=0 if(/Reaction Types/);
                next if(/^:/ or /^#/ or /^\s*$/);
                my @line_in_cfg=split ' ';
                if($line_in_cfg[1] eq "H2O") {
                    $line_in_cfg[1]="SPCF";
                }
                $map_from_resname_to_kernerl_number{$line_in_cfg[1]}=$line_in_cfg[2];
            }
        }
        close CFG;
        print %map_from_resname_to_kernerl_number if $debug==2;
        print "\n" if $debug==2;
        print "what are the resname , resid , segname?\n" if $debug;
    }
    elsif($input_count==2) {
        chomp;
        @resname_resid_segname=split ' ';
        my $count=0;
        for(@resname_resid_segname) {
            if($count % 3==2) {
                push @segnames,$_;
                print $_."\n" if $debug;
            }
            elsif($count % 3 ==1) {
                push @resids,$_;
                print $_." " if $debug;
            }
            else {
                push @resnames,$_;
                print $_." " if $debug;
            }
            $count++;
        }
	print "what name of evb.top?\n" if $debug;
    }
    elsif($input_count==3) {
	print if $debug;
	chomp;
	open TOP,">",$_;
        print "which pdb file?\n" if $debug;
    }
    elsif($input_count==4) {
        print if $debug;
        chomp;
        my $is_evb=0;$head_index=0;
	my $water_count=0;$hydronium_count=0;
        open PDB,"<",$_;
        
        my $atom_index_in_pdb=0;

        my $is_residue_start=0;
        my $is_residue_end=0;

        for(<PDB>) {
            next unless(/^ATOM/);
            my @line_in_pdb=split ' ';
            #my $atom_index_in_pdb=$line_in_pdb[1];
            $atom_index_in_pdb++;
	    last if $atom_index_in_pdb > $total_number_atoms;
	    if(/SPCF/) {
		print TOP sprintf "$atom_index_in_pdb\t%d\t%d\n",$map_from_resname_to_kernerl_number{"SPCF"},($water_count%3)+1;
		$water_count++;
		next;
	    }
	    elsif(/H3O/) {
		print TOP sprintf "$atom_index_in_pdb\t%d\t%d\n",$map_from_resname_to_kernerl_number{"H3O"},($hydronium_count%4)+1;
		$hydronium_count++;
		next;
	    }
            my $resname_in_pdb=$line_in_pdb[3];
            my $resid_in_pdb=$line_in_pdb[5];
            my $segname_in_pdb=$line_in_pdb[-2];
	   # print "$resname_in_pdb $resid_in_pdb $segname_in_pdb\n";
            for(my $i=0;$i<@resnames;$i++) { #THE FOLLOWING ALGO MAYBE BUGGY
                if($resnames[$i] eq $resname_in_pdb or $resnames[$i] eq $resname_in_pdb."P") {
                    if($resids[$i] eq $resid_in_pdb) {
                        if($segnames[$i] eq $segname_in_pdb) {
                            # make map_from_resname_to_kernerl_number work for GLUP ASPP ...:
                            $resname_in_pdb=$resnames[$i];
                            #IF HIT PERFECTLY
                            if($line_in_pdb[2] eq "N" and $is_residue_start==0) {
                                $is_residue_start=1;
                            }
                            elsif($line_in_pdb[2] eq "O" and $is_residue_start==1) {
                                $is_residue_start=0;
		                        print TOP sprintf "$atom_index_in_pdb\t%d\t%d\n",$map_from_resname_to_kernerl_number{$resname_in_pdb},$atom_index_in_pdb-$head_index+1;
                              $is_residue_end=1;
                            }
                            $is_evb=1;
			    last;
                        }
			else {$is_evb=0;}
                    }
		    else {$is_evb=0;}
                }
		else {$is_evb=0;}
	    }
        #if(($is_evb and $is_residue_start) or $line_in_pdb[2] eq "O") {
	    if(($is_evb and $is_residue_start)) { 
		$head_index=$atom_index_in_pdb if $head_index==0;
		print TOP sprintf "$atom_index_in_pdb\t%d\t%d\n",$map_from_resname_to_kernerl_number{$resname_in_pdb},$atom_index_in_pdb-$head_index+1;
	    }
	    else {
		$head_index=0;
      if(!$is_residue_end) {
		   print TOP "$atom_index_in_pdb\t0\t0\n";
      }
         $is_residue_end=0;
	    }
        }
        close PDB;
    }
    $input_count++;
    last if $input_count > 4;
}
close TOP;
