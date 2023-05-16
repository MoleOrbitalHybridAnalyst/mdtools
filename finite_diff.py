def finite_diff_grad(escanner, mol, dx=1e-5):
    import numpy as np
    coords0 = mol.atom_coords().copy()
    numgrad = np.zeros((mol.natm, 3))
    for i in range(mol.natm):
        for k in range(3):
            print(f"working on {i}-th atom {k}-th coordinate")
            coords = coords0.copy()
            coords[i][k] += dx / 2
            new_mol = mol.set_geom_(coords, unit='Bohr')
            Ep = escanner(new_mol)
            coords = coords0.copy()
            coords[i][k] -= dx / 2
            new_mol = mol.set_geom_(coords, unit='Bohr')
            Em = escanner(new_mol)
            numgrad[i][k] = (Ep - Em) / dx
    return numgrad
