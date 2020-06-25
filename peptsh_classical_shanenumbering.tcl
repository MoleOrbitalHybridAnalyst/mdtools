mol material AOShiny

mol showrep 0 0 off
mol addrep 0
mol modselect 1 0 protein
mol modstyle 1 0 NewCartoon
mol modmaterial 1 0 Transparent

mol addrep 0
mol modselect 2 0 name OW OH2 SOD CLA NA CL and z > 28 and z < 65 and within 5 of protein
mol modstyle 2 0 VDW
mol selupdate 2 0 1
mol modmaterial 1 0 Transparent

mol addrep 0
mol modselect 3 0 protein and resid 37 304 411 130 26 29 30 33
mol modstyle 3 0 Licorice

mol addrep 0
mol modselect 4 0 protein and resid 44 311 418 137 33 36 37 40
mol modstyle 4 0 Licorice
mol showrep 0 4 off

rotate x by -90
rotate y by 60
