#!/usr/bin/env python
import os
import numpy as np
import h5py
from pwscf_h5 import PwscfH5

if __name__ == '__main__':

    ref_fname = './ref/pwscf.pwscf.h5'
    ref = PwscfH5()
    ref.read(ref_fname)
    gvec = ref.get('gvectors')
    nkpt = int(ref.get('nkpt'))
    ref_kpts = np.zeros([nkpt,3])
    for ikpt in range(nkpt):
        rkpt = ref.val('electrons/kpoint_%d/reduced_k'%ikpt)
        ref_kpts[ikpt] = rkpt
    # end for
    ref_evals = ref.eigenvalues()

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

    # transfer orbital info.
    new.create_dataset('electrons/number_of_electrons',data=[1,1])
    new.create_dataset('electrons/number_of_kpoints',data=[8])
    new.create_dataset('electrons/number_of_spins',data=[1])
    new.create_dataset('electrons/psi_r_is_complex',data=[1])

    # transfer orbitals (electrons group)
    for ikpt in range(8): 
        kpt_path = 'electrons/kpoint_%d'%ikpt
        kgrp = new.create_group(kpt_path)
        kgrp.create_dataset('num_sym',data=[1])
        kgrp.create_dataset('symgroup',data=[1])
        kgrp.create_dataset('weight',data=[0.25])
        kgrp.create_dataset('reduced_k',data=ref_kpts[ikpt])
        if ikpt == 0:
            kgrp.create_dataset('gvectors',data=gvec)
            kgrp.create_dataset('number_of_gvectors',data=[len(gvec)])
        # end if 

        for ispin in range(1): # assume ispin==0
            spin_path = os.path.join(kpt_path,'spin_%d'%ispin)
            spgrp = new.create_group(spin_path)
            spgrp.create_dataset('number_of_states',data=[5]) # !!!! hard-code
            spgrp.create_dataset('eigenvalues',data=ref_evals[ikpt,ispin])
        
            for istate in range(4):
                loc = {'ikpt':ikpt,'ispin':ispin,'istate':istate}
                state_path = os.path.join(spin_path,'state_%d'%istate)

                # read, separate loop later
                psig_arr = ref.psig(**loc)

                # write
                psig_path = os.path.join(state_path,'psi_g')
                new.create_dataset(psig_path,data=psig_arr)
            # end for
    # end for

# end __main__
