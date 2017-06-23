#!/usr/bin/env python
from pwscf_h5 import PwscfH5
import numpy as np
import sys
sys.path.insert(0,'../../3-psir2hdf5/batch_convert')
from convert_psir_to_psig import isosurf

if __name__ == '__main__':

  vol = 78.7125735374 # simulation cell volume
  ikpt   = 0
  ispin  = 0
  istate = 0

  level_map = {0:0.01,1:0.015,2:0.01,3:0.015,4:0.012}
  level_map1 = {0:0.01,1:0.02,2:None,3:None,4:0.015}

  from pwscf_h5 import PwscfH5
  wf = PwscfH5()
  pyscf_h5 = '../1_workflow/pyscf2pwscf.h5'
  wf.read(pyscf_h5)
  wf_qe = PwscfH5()
  qe_h5 = '../../10-carbon-dimer/2_orbs/pwscf.pwscf.h5'
  # PBE0 calculation using exx_fraction=1.0
  #qe_h5 = '../../17-qe-hf/gamma/ecut20/scf/pwscf_output/pwscf.pwscf.h5'
  #rgrid_shape = np.array([16,16,16])*2+1
  wf_qe.read(qe_h5)
  for istate in range(4):
    psir = wf.get_psir_from_psig(ikpt,ispin,istate,rgrid_shape)
    val = np.absolute(psir)*vol/np.prod(rgrid_shape)

    psir1 = wf_qe.get_psir_from_psig(ikpt,ispin,istate,rgrid_shape)
    val1 = np.absolute(psir1)

    overlap = np.sum(val*val1)/np.sqrt(np.sum(val*val)*np.sum(val1*val1))
    print istate,overlap

    if overlap < 0.9999:
      import matplotlib.pyplot as plt
      fig = plt.figure()
      ax  = fig.add_subplot(121,projection='3d',aspect=1)
      ax1  = fig.add_subplot(122,projection='3d',aspect=1)

      # show orbitals
      ax.set_title('pyscf')
      ax1.set_title('qe')
      isosurf(ax,val)#,level_map[istate])
      isosurf(ax1,val1)#,level_map1[istate])

      # show carbon atoms
      pos = np.array([[0,0,0],[0.25,0.25,0.25]])*rgrid_shape # scale by grid size
      ax.plot(pos[:,0],pos[:,1],pos[:,2],'o',ms=20,c='y',alpha=0.6)
      ax1.plot(pos[:,0],pos[:,1],pos[:,2],'o',ms=20,c='y',alpha=0.6)
      #fig.savefig('diamond-carbon-2_orb%d.eps'%istate)
      plt.show()
    # end if
  # end for istate

# end __main__
