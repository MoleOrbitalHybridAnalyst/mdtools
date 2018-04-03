mol representation DynamicBonds 1.3 0.1 12.0
mol representation VDW 0.2 12.0
mol showrep 0 0 off
mol addrep 0
mol modselect 1 0 protein
mol modstyle 1 0 Licorice
mol addrep 0
mol modselect 2 0 not protein
mol modstyle 2 0 DynamicBonds
mol addrep 0
mol modselect 3 0 not protein
mol modstyle 3 0 VDW
mol addrep 0
mol modselect 4 0 serial 445
display projection orthographic
color Display Background white
color Labels Bonds black
source ~/share/vmdtools.tcl
# 899 961 962 963 966 967
