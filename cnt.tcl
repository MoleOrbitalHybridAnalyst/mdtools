mol representation DynamicBonds 1.4 0.1 12.0
mol representation VDW 0.2 12.0
mol showrep 0 0 off
mol addrep 0
mol modselect 1 0 carbon
mol modstyle 1 0 Lines
mol addrep 0
mol modselect 2 0 (oxygen or hydrogen) and within 5 of carbon
mol modstyle 2 0 DynamicBonds
mol selupdate 2 0 1
mol addrep 0
mol modselect 3 0 (oxygen or hydrogen) and within 5 of carbon
mol modstyle 3 0 VDW
mol selupdate 3 0 1
display projection orthographic
rotate y by 90
