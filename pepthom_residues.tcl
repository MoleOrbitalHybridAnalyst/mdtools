mol representation DynamicBonds 1.6 0.1 12.0
mol representation VDW 0.2 12.0
mol showrep 0 0 off
mol addrep 0
mol modselect 1 0 protein
mol modstyle 1 0 NewCartoon
mol showrep 0 1 off
mol addrep 0
mol modselect 2 0 resid 36 37 137 33 418 and protein
mol modstyle 2 0 Licorice
mol addrep 0
mol modselect 3 0 resid 44 311 79 40 41 and protein
mol modstyle 3 0 Licorice
##mol selupdate 3 0 1
##mol addrep 0
##mol modselect 4 0 resname H3O or same residue as (name OH2 and ((within 7 of (resid 799 and name OG1)) or (within 5 of (resid 767 and name OG)) or (within 5 of (resid 763 and name OH)) or (within 6 of (resid 908 and name OE1 OE2 OD1 OD2)) or (within 7 of (resid 794 and name NE1)) or (within 7 of (resid 944 and name ND1 NE2))))
##mol modstyle 4 0 VDW
##mol selupdate 4 0 1
###mol addrep 0
###mol modselect 4 0 resname H3O
###mol modstyle 4 0 CPK
##mol addrep 0
##mol modselect 5 0 resid 899 961 962 963 966 967 970 and protein
##mol modstyle 5 0 Licorice
##mol selupdate 5 0 1
##mol addrep 0
##mol modselect 6 0 same residue as (name OH2 and within 10 of index 14874)
##mol modstyle 6 0 DynamicBonds
##mol selupdate 6 0 1
display projection orthographic
##rotate x by 52.1348
##rotate z by 22.42915942
##rotate x by 90
# 899 961 962 963 966 967
