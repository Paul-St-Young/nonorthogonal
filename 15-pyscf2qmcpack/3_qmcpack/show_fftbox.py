#!/usr/bin/env python
def read_comp(fname):
  arr = []
  with open(fname,'r') as f:
    for line in f:
      tokens = line.strip('()\n').split(',')
      arr.append(map(float,tokens))
    # end for line
  # end with
  return np.array(arr,dtype=complex)
# end def

if __name__ == '__main__':
  vol = 78.7125735374
  level_map0= {0:0.01,1:0.015,2:0.01,3:0.015,4:0.012}
  level_map = {0:7,1:12,2:10,3:8}
  rgrid_shape = (9,9,9)

  from pwscf_h5 import PwscfH5
  wf = PwscfH5()
  wf.read('../1_workflow/pyscf2pwscf.h5')

  # get data
  import numpy as np
  for istate in range(4):
    rgrid = read_comp('fftbox%d.dat'%istate)
    psir  = (rgrid[:,0] + 1j*rgrid[:,1]).reshape(rgrid_shape)#*vol/np.prod(rgrid.shape)
    val  = np.absolute(psir)
    #print val.min(),val.max()
    psir0 = wf.get_psir_from_psig(0,0,istate,rgrid_shape=rgrid_shape)
    val0  = np.absolute(psir0)*vol/np.prod(rgrid_shape)

    overlap = np.sum(val0*val)/np.sqrt(np.sum(val*val)*np.sum(val0*val0))
    print istate,overlap

    # make plot
    import sys
    sys.path.insert(0,'../../6-verify_hdf5/3_check_no_ecut')
    from no_ecut_psig_to_moR import isosurf

    import matplotlib.pyplot as plt
    fig = plt.figure()
    ax  = fig.add_subplot(121,projection='3d',aspect=1)
    ax.set_title('pyscf hdf5')
    isosurf(ax,val0,level_map0[istate])
    ax1 = fig.add_subplot(122,projection='3d',aspect=1)
    ax1.set_title('QMCPACK fftbox')
    isosurf(ax1,val,level_map[istate])
    #plt.savefig('pyscf_qmcpack_orb%d.eps'%istate)
    plt.show()
  # end for istate

# end __main__
