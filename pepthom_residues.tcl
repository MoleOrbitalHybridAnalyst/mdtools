mol representation DynamicBonds 1.6 0.1 12.0
mol representation VDW 0.2 12.0
mol showrep 0 0 off
mol addrep 0
mol modselect 1 0 protein
mol modstyle 1 0 NewCartoon
mol showrep 0 1 off
mol addrep 0
#mol modselect 2 0 resid 44 311 and protein
mol modselect 2 0 resid 37 304 and protein
mol modstyle 2 0 Licorice
mol addrep 0
mol modselect 3 0 resid 26 29 30 126 and protein
mol modstyle 3 0 Licorice
mol addrep 0
mol modselect 4 0 (water and within 9 of (protein and residue 24)) or (water and within 4 of (protein and residue 126)) or (water and within 6 of (protein and residue 302))
mol modstyle 4 0 VDW
mol selupdate 4 0 1
mol addrep 0
mol modselect 5 0 (water and within 9 of (protein and residue 24)) or (water and within 4 of (protein and residue 126)) or resname H3O
mol modstyle 5 0 DynamicBonds
mol selupdate 5 0 1
display projection orthographic

# remove translation
proc remove_translation {} {
   set ref [atomselect 0 "backbone and residue 27 28 128 24 409 302 70 31 32" frame 0]
   set sel [atomselect 0 "backbone and residue 27 28 128 24 409 302 70 31 32"]
   set nf [molinfo get numframes]
   for {set i 0} {$i < $nf} {incr i} {
      $sel frame $i
      set trans [measure fit $sel $ref]
      [atomselect 0 all] move $trans
   }
}

# load the xtal pdb
mol new /project/gavoth/chhli/pot/pepthom_midway2/6exs.pdb
mol top 0
set ref [atomselect 0 "backbone and residue 27 28 128 24 409 302 70 31 32"]
set sel [atomselect 1 "backbone and residue 27 28 128 24 409 302 70 31 32"]
set trans [measure fit $sel $ref]
[atomselect 1 "all"] move $trans
mol addrep 1
mol modselect 1 1 resid 44 311 and protein
mol modstyle 1 1 Licorice
mol modmaterial 1 1 Transparent
mol addrep 1
mol modselect 2 1 resid 33 36 37 137 and protein
mol modstyle 2 1 Licorice
mol modmaterial 2 1 Transparent
mol addrep 1
mol modselect 3 1 water
mol modstyle 3 1 VDW
mol modmaterial 3 1 Transparent
mol showrep 1 0 off
mol showrep 1 3 off
