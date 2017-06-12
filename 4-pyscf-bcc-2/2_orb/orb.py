#!/usr/bin/env python
import os
import numpy as np
from pyscf.pbc import gto as pbcgto
import pyscf.pbc.tools.pyscf_ase as pyscf_ase

def build_cell(ase_atom,unit='B',ke=20.0,gsmax=None,basis='cc-pVDZ',
        pseudo=None,dimension=3):
    """ construct pyscf simulation cell from ase_atom """
    cell      = pbcgto.Cell()
    cell.unit = unit
    cell.atom = pyscf_ase.ase_atoms_to_pyscf(ase_atom)
    cell.a    = ase_atom.cell

    cell.basis  = basis
    cell.pseudo = pseudo
    cell.dimension = dimension

    if gsmax is not None: # different way to specify fft grid cutoff
        cell.gs = np.array([gsmax,gsmax,gsmax])
    else:
        cell.ke_cutoff = ke
    # end if

    cell.build()
    return cell
# end def build_cell

def run_dft(cell):
    from pyscf.pbc import dft as pbcdft
    sim = pbcdft.RKS(cell)
    sim.xc = 'lda,lda'
    sim.verbose = 3
    sim.scf()
    return sim
# end def run_dft

if __name__ == '__main__':

    from ase.build import bulk
    ase_atom = bulk('H','bcc',a=3.77945227,cubic=True)
    cell = build_cell(ase_atom,basis='cc-pVDZ')
    test = run_dft(cell)

    from pyscf.pbc.dft import gen_grid
    coords = gen_grid.gen_uniform_grids(test.cell)

    aoR = test._numint.eval_ao(cell,coords)
    moR = np.dot(aoR,test.mo_coeff)
    print test.mo_coeff.shape
    #moR = np.einsum('ri,ia->ra',aoR,test.mo_coeff)
    np.savetxt('moC.dat',test.mo_coeff)
    np.savetxt('moR.dat',moR)

    """ don't know how to pickle
    import pickle
    dft_obj = 'test.obj'
    if not os.path.isfile(dft_obj):
        test = run_dft(cell)
        with open(dft_obj,'w') as f:
            pickle.dump(test,f)
    else:
        with open(dft_obj,'r') as f:
            test = pickle.load(f)
    # end if
    """

# end __main__
