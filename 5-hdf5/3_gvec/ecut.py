#!/usr/bin/env python
import numpy as np
from pwscf_h5 import PwscfH5

def gvectors_within_cutoff(gs,ecut,recip_axes):
    """ e.g. gs=(4,4,4); ecut=20
     Inputs:
       ecut: kmax*kmax/2."""
    nx,ny,nz  = gs
    magsq_cut = 2*ecut

    # generate all gvector candidates
    from itertools import product
    all_int_gvecs = np.array([gvec for gvec in product(
      range(-nx,nx+1),range(-ny,ny+1),range(-nz,nz+1))])

    # make sure gvector magnitudes are within cutoff
    #  i.e. check g^2 < kmax^2 = 2*ecut
    all_gvecs = np.dot(all_int_gvecs,recip_axes)
    all_mag2  = np.sum(all_gvecs*all_gvecs,axis=1)
    sel = (all_mag2 < magsq_cut)
    return all_int_gvecs[sel]
# end def gvectors_within_cutoff

if __name__ == '__main__':

    ref = PwscfH5()
    ref.read('../1_transfer_old/ref/pwscf.pwscf.h5')

    # get axes
    axes = ref.get('axes')
    recip_axes = 2*np.pi*np.linalg.inv(axes)

    #gvecs = gvectors_within_cutoff((2,2,2),20,recip_axes)

    # get bcc-2 cell
    import sys
    sys.path.insert(0,'../2_system_from_cell')
    from system_from_cell import bcc2
    alat = 3.77945227
    cell = bcc2(alat,ke=20)
    assert np.allclose(axes,cell.a)
    vol = np.dot(np.cross(axes[0],axes[1]),axes[2])

    # put AO on numerical grid
    from pyscf.pbc.dft import gen_grid,numint
    coords = gen_grid.gen_uniform_grids(cell)
    aoR    = numint.eval_ao(cell,coords)

    # double-check that I remember how to use np.fft
    #  transform the first atomic orbital back and forth
    aoR0 = aoR[:,0].reshape(9,9,9)
    #aoG0 = np.fft.fftn(np.fft.fftshift(aoR0))/np.prod(aoR0.shape)*vol
    #aoR1 = np.fft.fftshift(np.fft.ifftn(aoG0))*np.prod(aoR0.shape)/vol
    aoG0 = np.fft.fftn(aoR0)/np.prod(aoR0.shape)*vol
    aoR1 = np.fft.ifftn(aoG0)*np.prod(aoR0.shape)/vol

    if not np.allclose(aoR0,aoR1):
        sys.path.insert(0,'../../3-psir2hdf5/batch_convert')
        from convert_psir_to_psig import isosurf
        import matplotlib.pyplot as plt
        fig = plt.figure()
        ax  = fig.add_subplot(121,projection='3d',aspect=1)
        ax1 = fig.add_subplot(122,projection='3d',aspect=1)
        isosurf(ax,aoR0,0.2)
        isosurf(ax1,aoR1,0.2)
        plt.show()

        raise AssertionError('fft failed!')
    # end if
    # note: forward transform normalization is *vol/ngrid
    #      backward transform normalization is *ngrid/vol
    #  no fftshift is needed for transform
    fft_norm = vol/np.prod(aoR0.shape)

    # guess the gvectors
    nx,ny,nz  = (4,4,4)
    from itertools import product
    all_int_gvecs = np.array([gvec for gvec in product(
      range(-nx,nx+1),range(-ny,ny+1),range(-nz,nz+1))])
    all_gvecs = np.dot(all_int_gvecs,recip_axes)
    kgrid = np.zeros(aoG0.shape,dtype=complex)
    #rgrid = np.zeros(aoR0.shape,dtype=complex)
    for igvec in range(len(all_gvecs)):
      int_gvec = all_int_gvecs[igvec]
      gvec = all_gvecs[igvec]
      kgrid[tuple(int_gvec)] = aoG0[tuple(int_gvec)]
      #rvec = np.dot(int_gvec,axes)/aoR0.shape
      #rgrid[tuple(int_gvec)] += np.exp(1j*np.dot(gvec,rvec)*2*np.pi)*aoG0[tuple(int_gvec)]
    # end for igvec
    test = np.fft.ifftn(kgrid)/fft_norm
    assert np.allclose(test,aoR0)

    """
    print axes

    # get gvectors
    gvec = ref.get('gvectors')
    abs_gvec = np.dot(gvec,recip_axes)
    print gvec

    # check maximum kinetic
    gmag = np.linalg.norm(abs_gvec,axis=1)
    gmax = gmag.max()
    print (gmax*gmax).sum()
    """
