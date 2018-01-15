mol showrep 0 0 off
mol addrep 0
mol modselect 1 0 same residue as (type OH2 OW OT and  within 7 of (resid 425 322 and name OE1 OE2 OD1 OD2))
mol modstyle 1 0 bonds
mol addrep 0
mol modselect 2 0 same residue as (type OH2 OW OT and within 4 of (resid 37 and name NH1 NH2) )
mol modstyle 2 0 bonds
mol addrep 0
mol modselect 3 0 resid 425 322 324 37 and protein
mol modstyle 3 0 bonds
mol addrep 0
mol modselect 4 0 same residue as (type OH2 OW OT and  within 4 of (resid 324 and name NZ))
mol modstyle 4 0 bonds
mol addrep 0
mol modselect 5 0 segname W3
mol modstyle 5 0 bonds
mol addrep 0
mol modselect 6 0 type CLA
mol modstyle 6 0 vdw
#mol addrep 0
#mol modselect 5 0 waters and x<51.5 and x>37.5 and y>38 and y<51.5 and z<58 and z>38.5
#mol modstyle 5 0 bonds
display projection orthographic
rotate x by 90
rotate z by 180
