#!/usr/bin/env python
import os
import numpy as np
import pandas as pd
import h5py
from pwscf_h5 import PwscfH5

if __name__ == '__main__':

    ref_fname = './ref/pwscf.pwscf.h5'
    ref = PwscfH5()
    ref.read(ref_fname)
    gvec = ref.get('gvectors')
    nkpt0 = ref.get('nkpt')[0]
    ref_evals = ref.eigenvalues()
    nkpt,nspin,nstate = ref_evals.shape
    assert nkpt == nkpt0
    ref_kpts = np.zeros([nkpt,3])
    for ikpt in range(nkpt):
        rkpt = ref.val('electrons/kpoint_%d/reduced_k'%ikpt)
        ref_kpts[ikpt] = rkpt
    # end for

    new_fname = 'from_ref.h5'
    new = h5py.File(new_fname,'w')

    # transfer metadata
    # ====
    # transfer version info.
    new.create_dataset('application/code',data=['espresso'])
    new.create_dataset('application/version',data=[4,1,4])
    new.create_dataset('format',data=['ES-HDF'])
    new.create_dataset('version',data=[2,1,0])

    # transfer lattice (supercell group)
    alat = 3.77945227
    axes = alat*np.eye(3)
    new.create_dataset('supercell/primitive_vectors',data=axes)

    # transfer particle set (atoms group)
    atoms = ['H','H']
    new.create_dataset('atoms/number_of_atoms',data=[len(atoms)])

    pos   = np.dot(np.array([[0,0,0],[0.5,0.5,0.5]]), axes)
    assert len(atoms)==len(pos)
    new.create_dataset('atoms/positions',data=pos)

    species = np.unique(atoms)
    new.create_dataset('atoms/number_of_species',data=[len(species)])

    species_map = {'H':0,'He':1,'Li':2,'Be':3}
    species_ids = [species_map[atom] for atom in atoms]
    new.create_dataset('atoms/species_ids',data=species_ids)

    # describe existing species (!!!! just do H for now)
    sp0 = new.create_group('atoms/species_0')
    sp0.create_dataset('name',data=['H'])
    sp0.create_dataset('atomic_number',data=[1])
    sp0.create_dataset('valence_charge',data=[1])

    # transfer wavefunction data
    # ====
    ref.create_electrons_group(new,gvec,ref.eigensystem(),[1,1])

# end __main__
