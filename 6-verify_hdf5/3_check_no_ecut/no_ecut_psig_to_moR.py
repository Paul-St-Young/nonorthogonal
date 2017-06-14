#!/usr/bin/env python
import numpy as np
from pwscf_h5 import PwscfH5
import sys
sys.path.insert(0,'../../3-psir2hdf5/batch_convert/')
from convert_psir_to_psig import isosurf
sys.path.insert(0,'../1_psig_to_moR/')
from psig_to_moR import simulate_bcc2, psig_to_psir, get_psir, get_pyscf_psir

def check_orb(pwscf_h5_fname,ikpt=0,ispin=0,istate=0,thres=5):
    """ check  """
    psir0 = get_pyscf_psir(kmf,cell,**loc)
    val0  = np.absolute(psir0)

    # compare to pwscf.pwscf.h5
    wf = PwscfH5()
    wf.read(pwscf_h5_fname)
    rkpt_path = 'electrons/kpoint_%d/reduced_k'%loc['ikpt']
    rkpt  = wf.fp[rkpt_path].value
    assert np.allclose(rkpts0[loc['ikpt']],rkpt)
    psir  = get_psir(fname,**loc)
    val  = np.absolute(psir)

    # see if orbitals are the same
    percent_difference = np.linalg.norm(val0-val)/np.linalg.norm(val0)*100.
    if percent_difference > thres:
      #print val0.min(),val0.max()
      #print val.min(),val.max()
      print 'ikpt=%d istate=%d'%(ikpt,istate), percent_difference,'% different'
      import matplotlib.pyplot as plt
      fig = plt.figure()
      ax0 = fig.add_subplot(121,projection='3d',aspect=1)
      ax  = fig.add_subplot(122,projection='3d',aspect=1)
      isosurf(ax0,val0)
      isosurf(ax,val)
      plt.show()
    # end if
# end def check_orb

if __name__ == '__main__':

    # run DFT
    kgrid = [2,2,2]
    cell,kmf = simulate_bcc2(
      alat=3.77945227, basis='cc-pVDZ',chkfile='pvdz.h5',ke=20.,kgrid=kgrid)
    abs_kpts = cell.make_kpts(kgrid)
    rkpts0   = cell.get_scaled_kpts(abs_kpts)

    # define orbital to compare
    fname ='../2_no_ecut_hdf5/pwscf.pwscf.h5'
    for ikpt in range(8):
      for istate in range(5):
        print ikpt,istate
        loc = {'ikpt':ikpt,'ispin':0,'istate':istate}
        check_orb(fname,thres=0.1,**loc)
      # end for
    # end for

# end __main__
