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

  # get data
  import numpy as np
  #rgrid = np.genfromtxt('fftbox.dat',delimiter=',',dtype=np.complex128)
  for istate in range(4):
    rgrid = read_comp('fftbox%d.dat'%istate)
    psir  = (rgrid[:,0] + 1j*rgrid[:,1]).reshape(9,9,9)
    val  = np.absolute(psir)
    print val.min(),val.max()

    # make plot
    import sys
    sys.path.insert(0,'../../6-verify_hdf5/3_check_no_ecut')
    from no_ecut_psig_to_moR import isosurf

    import matplotlib.pyplot as plt
    fig = plt.figure()
    ax  = fig.add_subplot(111,projection='3d',aspect=1)
    isosurf(ax,val)
    plt.show()

# end __main__
