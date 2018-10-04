mol representation DynamicBonds 1.6 0.1 12.0
mol representation VDW 0.2 12.0
mol showrep 0 0 off
mol addrep 0
mol modselect 1 0 protein
mol modstyle 1 0 NewCartoon
mol showrep 0 1 off
mol addrep 0
mol modselect 2 0 residue 27 28 128 24 409 and protein
mol modstyle 2 0 Licorice
mol addrep 0
mol modselect 3 0 residue 35 302 70 31 32 and protein
mol modstyle 3 0 Licorice
mol addrep 0
mol modselect 4 0 (water and within 9 of (protein and residue 24)) or (water and within 4 of (protein and residue 132))
mol modstyle 4 0 VDW
mol selupdate 4 0 1
mol addrep 0
mol modselect 5 0 (water and within 9 of (protein and residue 24)) or (water and within 4 of (protein and residue 132))
mol modstyle 5 0 DynamicBonds
mol selupdate 5 0 1
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
