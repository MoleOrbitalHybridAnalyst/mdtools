mol representation DynamicBonds 1.4 0.1 12.0
mol representation VDW 0.2 12.0
mol showrep 0 0 off
mol addrep 0
mol modselect 1 0 carbon
mol modstyle 1 0 Lines
mol addrep 0
mol modselect 2 0 oxygen or hydrogen within 5 of carbon
mol modstyle 2 0 DynamicBonds
mol addrep 0
mol modselect 3 0 oxygen or hydrogen within 5 of carbon
mol modstyle 3 0 VDW
display projection orthographic
rotate y by 90
