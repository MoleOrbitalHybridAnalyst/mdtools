mol modselect 0 0 protein
mol modstyle 0  0  NewCartoon
mol showrep 0 0 off
mol addrep 0
mol modselect 1 0 resid 67 322 425 460 and protein
mol modstyle 1 0 bonds
mol addrep 0
mol modselect 2 0 same residue as type OW OH2 OTH and within 6 of (resid 67 and protein)
mol modstyle 2 0 cpk

set cec [atomselect 0 "index 29000"]
$cec moveto {40.054857 39.764586 36.782492}

display projection orthographic
rotate x by 90
