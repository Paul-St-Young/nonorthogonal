#!/usr/bin/env python
import os
import numpy as np
import sys
sys.path.insert(0,'../2_orb/')
from orb import build_cell

def run_kdft(cell,nmp=(2,2,2),gamma=True):
    import ase
    from pyscf.pbc import dft as pbcdft

    # define kpoints in reciprocal lattice units
    scaled_kpts = np.array([
        [0,0,0],
        [0,0,-0.5],
        [0,-0.5,0],
        [0,-0.5,-0.5],
        [-0.5,0,0],
        [-0.5,0,-0.5],
        [-0.5,-0.5,0],
        [-0.5,-0.5,-0.5]
    ])

    # get kpoints in momentum units
    inv_lat  = np.linalg.inv(cell.lattice_vectors())
    abs_kpts = np.dot(scaled_kpts,inv_lat)
    #np.savetxt('abs_kpts.dat',abs_kpts)
    #abs_kpts = cell.get_abs_kpts(scaled_kpts)

    # make simulation
    kmf = pbcdft.KRKS(cell, abs_kpts)
    # kpoint-average mean-filed (kmf) calculation

    kmf.xc = 'lda,lda'
    kmf.verbose = 3
    kmf.scf()
    return kmf
# end def run_dft

if __name__ == '__main__':

    from ase.build import bulk
    ase_atom = bulk('H','bcc',a=3.77945227,cubic=True)
    cell = build_cell(ase_atom,basis='cc-pVDZ',ke=20)
    #test = run_kdft(cell,nmp=[4,4,4],gamma=True)

    #from pyscf.pbc.dft import gen_grid
    #coords = gen_grid.gen_uniform_grids(test.cell)
    #aoR = test._numint.eval_ao(cell,coords)[0]

    #np.savetxt('aoR.dat',aoR)

    """
    coeff = test.mo_coeff
    print coeff.shape
    #print coeff[0,:,:]/coeff[1,:,:]
    #assert np.allclose(coeff[0,:,:],coeff[1,:,:])
    data = np.zeros([coeff.shape[0],aoR.shape[0],coeff.shape[1]],dtype=complex)
    for ikpt in range(8):
        data[ikpt,:,:] = np.dot(aoR,coeff[ikpt,:,:])
    # end for ikpt
    #moR = np.dot(aoR,test.mo_coeff.transpose([1,2,0]))
    #data = moR.transpose([2,0,1]).reshape(8,729*10)
    #np.savetxt('moR.dat',data.real.reshape(8,729*10))
    import pickle
    with open('data.p','w') as f:
        pickle.dump(data,f)
    """

# end __main__
