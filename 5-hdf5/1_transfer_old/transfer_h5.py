#!/usr/bin/env python
import os
import h5py
from pwscf_h5 import PwscfH5

if __name__ == '__main__':

    ref_fname = './ref/pwscf.pwscf.h5'
    ref = PwscfH5()
    ref.read(ref_fname)
    gvec = ref.get('gvectors')

    new_fname = 'from_ref.h5'
    new = h5py.File(new_fname,'w')

    kpt0 = new.create_group('electrons/kpoint_0')
    kpt0.create_dataset('gvectors',data=gvec)
    
    for ikpt in range(8): # assume ispin==0
        for istate in range(4):
            loc = {'ikpt':ikpt,'ispin':0,'istate':istate}
            psig_arr = ref.psig(**loc)

            psig_path = 'electrons/kpoint_{kpt:d}/spin_0/state_{state:d}/psi_g'.format(kpt=ikpt,state=istate)
            new.create_dataset(psig_path,data=psig_arr)
        # end for
    # end for
