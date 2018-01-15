mol showrep 0 0 off
mol addrep 0
mol modselect 1 0 protein
mol modstyle 1 0 NewCartoon
mol showrep 0 1 off
mol addrep 0
mol modselect 2 0 resid 799 798 794 908 905 944 and protein
mol modstyle 2 0 Licorice
mol addrep 0
mol modselect 3 0 waters and name OH2 and z <=123 and z>=107
mol modstyle 3 0 VDW
display projection orthographic
rotate x by -90

